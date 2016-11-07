import pysam
import numpy as np
import ctypes
from math import ceil

from multiprocessing import Process, Queue, Array, Value
import sys
import multiprocessing, logging

mpl = multiprocessing.log_to_stderr()
mpl.setLevel(logging.INFO)

def partition_snps(region, n_parts, start_1pos, end_1pos):
    snps_in_window = sum(region) / float(n_parts-1)
    curr_window = 0
    window_positions = [start_1pos]
    for curr_0pos, flag in enumerate(region[start_1pos-1 : end_1pos]):
        curr_window += flag
        if curr_window >= snps_in_window:
            window_positions.append(curr_0pos+1) #TODO HELP
            curr_window = 0
    return window_positions

#TODO SENTINEL SYMBOLS BEFORE AND AFTER A READ
#TODO What happens if we traverse backwards...?
#TODO Single SNP reads could use a pairwise observation with themselves? (A, A, i, i)
def load_from_bam(h, bam, target_contig, start_pos, end_pos, vcf_handler, use_end_sentinels=False, n_threads=1):
    """
    Load variants observed in a :py:class:`pysam.AlignmentFile` to
    an instance of :py:class:`hansel.hansel.Hansel`.

    Parameters
    ----------
    hansel : :py:class:`hansel.hansel.Hansel`
        An initialised instance of the `Hansel` data structure.

    bam : :py:class:`pysam.AlignmentFile`
        A BAM alignment.

    target_contig : str
        The name of the contig for which to recover haplotypes.

    start_pos : int
        The 1-indexed genomic position from which to begin considering variants.

    end_pos : int
        The 1-indexed genomic position at which to stop considering variants.

    vcf_handler : dict{str, any}
        Variant metadata, as provided by :py:func:`gretel.gretel.process_vcf`.

    use_end_sentinels : boolean, optional(default=False)
        Whether or not to append an additional pairwise observation between
        the final variant on a read towards a sentinel.

        .. note:: Experimental
          This feature is for testing purposes, currently it is recommended
          that the flag be left at the default of `False`. However, some
          data sets report minor performance improvements for some haplotypes
          when set to `True`. Currently there is no applicable re-weighting
          scheme for reducing the observations that end at sentinels.
          This flag may be removed at any time without warning.

    n_threads : int, optional(default=1)
        Number of threads to spawn for reading the BAM

    Returns
    -------
    Metadata : dict{str, any}
        A dictionary of metadata that may come in useful later.
        Primarily used to return a list of integers describing the number of
        variants covered by each read in the provided alignment BAM.

    Raises
    ------
    Exception
        Aborts if an unsupported CIGAR operation is encountered.
        Currently supports Match, Insert, Delete, Soft Clip.
    Exception
        Aborts if a SNP site is not reached when parsing the CIGAR for a read that
        should overlap that position.
    """

    meta = {}


    hansel = np.frombuffer(Array(ctypes.c_float, 6 * 6 * (vcf_handler["N"]+2) * (vcf_handler["N"]+2), lock=False), dtype=ctypes.c_float)
    hansel = hansel.reshape(6, 6, vcf_handler["N"]+2, vcf_handler["N"]+2)
    hansel.fill(0.0)

    import random
    def progress_worker(progress_q, n_workers, slices, total_snps):
        worker_pos = []
        worker_done = []
        for _ in range(0, n_workers):
            worker_pos.append(0)
            worker_done.append(0)

        while sum(worker_done) < n_workers:
            work_block = progress_q.get()
            worker_pos[work_block["worker_i"]] = work_block["pos"]
            if not work_block["pos"]:
                worker_done[work_block["worker_i"]] = 1

                slices.value += work_block["slices"]
                total_snps.value += work_block["covered_snps"]
                sys.stderr.write("%s\n" % ([ worker_pos[i] if status != 1 else None for (i, status) in enumerate(worker_done)]))
            if random.random() < 0.1:
                sys.stderr.write("%s\n" % ([ worker_pos[i] if status != 1 else None for (i, status) in enumerate(worker_done)]))
        return (slices, total_snps)

    def bam_worker(bam_q, progress_q, worker_i):
        def __symbol_num(symbol):
            symbols = ['A', 'C', 'G', 'T', 'N', '_']
            #TODO Catch potential IndexError
            #TODO Generic mechanism for casing (considering non-alphabetically named states, too...)
            return symbols.index(symbol)

        symbols = ['A', 'C', 'G', 'T', 'N', '_']
        unsymbols = ['_', 'N']
        worker = worker_i

        slices = 0
        covered_snps = 0

        while True:
            work_block = bam_q.get()
            if work_block is None:
                progress_q.put({
                    "pos": None,
                    "worker_i": worker_i,
                    "slices": slices,
                    "covered_snps": covered_snps,
                })
                break

            for read in bam.fetch(target_contig, start=work_block["start"]-1, end=work_block["end"], multiple_iterators=True):

                START_POS_OFFSET = 0

                if read.is_duplicate or read.is_secondary:
                    continue

                LEFTMOST_1pos = read.reference_start + 1 # Convert 0-based reference_start to 1-based position (to match region array and 1-based VCF)
                RIGHTMOST_1pos = read.reference_end #ofc this is 1-indexed instead of 0

                # Special case: Consider reads that begin before the start_pos, but overlap the 0th block
                if work_block["i"] == 0:
                    if LEFTMOST_1pos < start_pos:
                        # Read starts before the start_pos
                        if read.reference_start + 1 + read.query_alignment_length < start_pos:
                            # Read ends before the start_pos
                            continue

                        LEFTMOST_1pos = start_pos
                        START_POS_OFFSET = (start_pos - (read.reference_start + 1))

                else:
                    # This read begins before the start of the current (non-0) block
                    # and will have already been covered by the block that preceded it
                    if LEFTMOST_1pos < work_block["start"]:
                        continue

                # If the current read begins after the region of interest, stop parsing the sorted BAM
                #if LEFTMOST_1pos > end_pos:
                #    break

                # Read ends after the end_pos of interest, so clip it
                if RIGHTMOST_1pos > end_pos:
                    RIGHTMOST_1pos = end_pos

                # Check if the read actually covers any SNPs
                support_len = np.sum(vcf_handler["region"][LEFTMOST_1pos : RIGHTMOST_1pos + 1])

                # Ignore reads without evidence
                if support_len == 0:
                    continue

                slices += 1

                rank = np.sum(vcf_handler["region"][1 : LEFTMOST_1pos])
                support_seq = []
                for i in range(0, support_len):
                    snp_rev = vcf_handler["snp_rev"][rank + i]

                    snp_pos_on_read = snp_rev - LEFTMOST_1pos + START_POS_OFFSET

                    aligned_residues = [x for x in read.get_aligned_pairs(with_seq=True) if x[1] is not None] # Filter out SOFTCLIP and INS
                    snp_pos_on_aligned_read = aligned_residues[snp_pos_on_read][0]

                    try:
                        support_seq.append(read.query_sequence[snp_pos_on_aligned_read])
                    except TypeError:
                        sys.stderr.write("NoneType (DEL) found at SNP site on read '%s', reference position %d\n" % (read.qname, snp_rev))
                        support_seq.append("_")

                support_seq = "".join(support_seq)
                covered_snps += len(support_seq.replace("N", "").replace("-", ""))

                progress_q.put({"pos": read.reference_start + 1, "worker_i": worker_i})

                # For each position in the supporting sequence (that is, each covered SNP)
                for i in range(0, support_len):
                    snp_a = support_seq[i]

                    #if support_len == 1:
                    #    if rank == 0:
                    #        hansel.add_observation('_', snp_a, 0, 1)
                    #        hansel.add_observation(snp_a, '_', 1, 2)
                    #    else:
                    #        hansel.add_observation(snp_a, '_', rank+1, rank+2)


                    # For each position in the supporting sequence following i
                    for j in range(i+1, support_len):
                        snp_b = support_seq[j]

                        # Ignore observations who are from an invalid transition
                        if snp_a in unsymbols:
                            continue

                        # Sentinel->A
                        if i==0 and j==1 and rank==0:
                            # If this is the first position in the support (support_pos == 0)
                            # and rank > 0 (that is, this is not the first SNP)
                            # and SNPs a, b are adjacent
                            hansel[__symbol_num('_'), __symbol_num(snp_a), 0, 1] += 1
                            hansel[__symbol_num(snp_a), __symbol_num(snp_b), 1, 2] += 1

                        # B->Sentinel
                        elif (j+rank+1) == vcf_handler["N"] and abs(i-j)==1:
                            # Last observation (abs(i-j)==1),
                            # that ends on the final SNP (j+rank+1 == N)
                            hansel[__symbol_num(snp_a), __symbol_num(snp_b), vcf_handler["N"]-1, vcf_handler["N"]] += 1
                            hansel[__symbol_num(snp_b), __symbol_num('_'), vcf_handler["N"], vcf_handler["N"]+1] += 1

                        # A regular observation (A->B)
                        else:
                            hansel[__symbol_num(snp_a), __symbol_num(snp_b), i+rank+1, j+rank+1] += 1

                            if use_end_sentinels:
                                if j==(support_len-1) and abs(i-j)==1:
                                    # The last SNP on a read, needs a sentinel afterward
                                    hansel[__symbol_num(snp_b), __symbol_num('_'), j+rank+1, j+rank+2] += 1

    bam_queue = Queue()
    progress_queue = Queue()

    # Queue the wokers
    # TODO Evenly divide, but in future, consider the distn
    window_l = int((end_pos - start_pos) / float(n_threads))
    for window_i, window_pos in enumerate(range(start_pos, end_pos+1, window_l)):
        bam_queue.put({
            "start": window_pos,
            "end": window_pos + window_l,
            "i": window_i,
        })

    processes = []
    for _ in range(n_threads):
        p = Process(target=bam_worker,
                    args=(bam_queue, progress_queue, _))
        processes.append(p)

    # ...and a progress process
    n_reads = Value('i', 0)
    total_covered_snps = Value('i', 0)
    p = Process(target=progress_worker,
                args=(progress_queue, n_threads, n_reads, total_covered_snps))
    processes.append(p)

    for p in processes:
        p.start()

    # Add sentinels
    for _ in range(n_threads):
        bam_queue.put(None)

    # Wait for processes to complete work
    for p in processes:
        p.join()

    meta["hansel"] = hansel
    meta["L"] = int(ceil(float(total_covered_snps.value)/n_reads.value))
    return meta

def load_fasta(fa_path):
    """
    A convenient wrapper function for constructing a :py:class:`pysam.FastaFile`

    Parameters
    ----------
    fa_path : str
        Path to FASTA

    Returns
    -------

    FASTA File Interface : :py:class:`pysam.FastaFile`
    """
    return pysam.FastaFile(fa_path)

