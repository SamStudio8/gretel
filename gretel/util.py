import pysam
import numpy as np
from math import ceil
from hansel import Hansel
import vcf

from multiprocessing import Process, Queue, Value
import sys

def get_ref_len_from_bam(bam_path, target_contig):
    """
    Fetch the length of a given reference sequence from a :py:class:`pysam.AlignmentFile`.

    Parameters
    ----------
    bam_path : str
        Path to the BAM alignment

    target_contig : str
        The name of the contig for which to recover haplotypes.

    Returns
    -------
    end_pos : int
        The 1-indexed genomic position at which to stop considering variants.
    """
    bam = pysam.AlignmentFile(bam_path)
    end = bam.lengths[bam.get_tid(target_contig)]
    bam.close()

    return end

def load_from_bam(bam_path, target_contig, start_pos, end_pos, vcf_handler, use_end_sentinels=False, n_threads=1, debug_reads=False, debug_pos=False):
    """
    Load variants observed in a :py:class:`pysam.AlignmentFile` to
    an instance of :py:class:`hansel.hansel.Hansel`.

    Parameters
    ----------
    bam_path : str
        Path to the BAM alignment

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
          when set to `True`.
          This flag may be removed at any time without warning.

    n_threads : int, optional(default=1)
        Number of threads to spawn for reading the BAM

    Returns
    -------
    Hansel : :py:class:`hansel.hansel.Hansel`
    """

    hansel = Hansel.init_matrix(['A', 'C', 'G', 'T', 'N', "-", "_"], ['N', "_"], vcf_handler["N"])

    if not debug_reads:
        debug_reads = set([])
    if not debug_pos:
        debug_pos = set([])

    import random
    def progress_worker(progress_q, n_workers, slices, total_snps, crumbs):
        worker_pos = []
        worker_done = []
        for _ in range(0, n_workers):
            worker_pos.append(0)
            worker_done.append(0)

        while sum(worker_done) < n_workers:
            work_block = progress_q.get()
            worker_pos[work_block["worker_i"]] = work_block["pos"]
            if work_block["pos"] is None:
                worker_done[work_block["worker_i"]] = 1

                crumbs.value += work_block["crumbs"]
                slices.value += work_block["slices"]
                total_snps.value += work_block["covered_snps"]
                sys.stderr.write("%s\n" % ([ worker_pos[i] if status != 1 else None for (i, status) in enumerate(worker_done)]))
            if random.random() < 0.1:
                sys.stderr.write("%s\n" % ([ worker_pos[i] if status != 1 else None for (i, status) in enumerate(worker_done)]))
        return (slices, total_snps, crumbs)

    def bam_worker(bam_q, progress_q, worker_i):

        worker = worker_i

        slices = 0
        crumbs = 0
        covered_snps = 0

        bam = pysam.AlignmentFile(bam_path)

        while True:
            work_block = bam_q.get()
            if work_block is None:
                progress_q.put({
                    "pos": None,
                    "worker_i": worker_i,
                    "slices": slices,
                    "crumbs": crumbs,
                    "covered_snps": covered_snps,
                })
                break

            reads = {}
            for p_col in bam.pileup(reference=target_contig, start=work_block["start"]-1, end=work_block["end"], ignore_overlaps=False, min_base_quality=0):

                if p_col.reference_pos + 1 > end_pos:
                    # Ignore positions beyond the end_pos
                    break

                if vcf_handler["region"][p_col.reference_pos+1] != 1:
                    continue

                for p_read in p_col.pileups:

                    curr_read_1or2 = 0
                    if p_read.alignment.is_paired:
                        if p_read.alignment.is_read1:
                            curr_read_1or2 = 1
                        elif p_read.alignment.is_read2:
                            curr_read_1or2 = 2
                        else:
                            #TODO Probably indicative of bad data
                            pass


                    curr_read_name = "%s_%d" % (p_read.alignment.query_name, curr_read_1or2)

                    LEFTMOST_1pos = p_read.alignment.reference_start + 1 # Convert 0-based reference_start to 1-based position (to match region array and 1-based VCF)

                    # Special case: Consider reads that begin before the start_pos, but overlap the 0th block
                    if work_block["i"] == 0:
                        if LEFTMOST_1pos < start_pos:
                            # Read starts before the start_pos
                            if p_read.alignment.reference_start + 1 + p_read.alignment.query_alignment_length < start_pos:
                                # Read ends before the start_pos
                                continue
                            LEFTMOST_1pos = start_pos
                            #continue
                    else:
                        # This read begins before the start of the current (non-0) block
                        # and will have already been covered by the block that preceded it
                        if LEFTMOST_1pos < work_block["start"]:
                            continue

                    sequence = None
                    qual = None
                    if p_read.is_del:
                        # TODO Not sure about how to estimate quality of deletion?
                        sequence = "-" * (abs(p_read.indel) + 1)
                        qual = p_read.alignment.query_qualities[p_read.query_position_or_next] * (abs(p_read.indel) + 1)
                    elif p_read.indel > 0:
                        # p_read.indel peeks to next CIGAR and determines whether the base FOLLOWING this one is an insertion or not
                        sequence = p_read.alignment.query_sequence[p_read.query_position : p_read.query_position + p_read.indel + 1]
                        qual = p_read.alignment.query_qualities[p_read.query_position : p_read.query_position + p_read.indel + 1]
                    else:
                        sequence = p_read.alignment.query_sequence[p_read.query_position]
                        qual = p_read.alignment.query_qualities[p_read.query_position]

                    if not sequence:
                        print("[WARN] Sequence data seems to not be correctly salvaged from read %s" % p_read.alignment.query_name)
                        continue

                    if curr_read_name not in reads:
                        reads[curr_read_name] = {
                            "rank": np.sum(vcf_handler["region"][1 : LEFTMOST_1pos]),
                            "seq": [],
                            "quals": [],
                            "refs_1pos": [],
                            "read_variants_0pos": [],
                        }
                    reads[curr_read_name]["seq"].append(sequence)
                    reads[curr_read_name]["quals"].append(qual)
                    reads[curr_read_name]["refs_1pos"].append(p_col.reference_pos+1)
                    reads[curr_read_name]["read_variants_0pos"].append(p_read.query_position)


            for dread in debug_reads:
                for read_type in [0,1,2]:
                    r = reads.get((dread+"_"+str(read_type)),None)
                    if r:
                        for snp_i, ref_pos in enumerate(r["refs_1pos"]):
                            print (dread, read_type, ref_pos, r["seq"][snp_i])
                        print ("RANK", dread, read_type, r["rank"])

            if debug_pos:
                for read in reads:
                    for d_pos in set(reads[read]["refs_1pos"]) & debug_pos:
                        i = reads[read]["refs_1pos"].index(d_pos)
                        print (read, d_pos, reads[read]["seq"][i])


            num_reads = len(reads)
            for qi, qname in enumerate(reads):
                progress_q.put({"pos": num_reads-(qi+1), "worker_i": worker_i})

                if not len(reads[qname]["seq"]) > 1:
                    # Ignore reads without evidence
                    continue
                slices += 1

                rank = reads[qname]["rank"]
                support_len = len(reads[qname]["seq"])

                support_seq = "".join([b[0] for b in reads[qname]["seq"]]) # b[0] has the affect of capturing the base before any insertion
                covered_snps += len(support_seq.replace("N", "").replace("_", ""))

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
                        if snp_a in ['_', 'N']:
                            continue

                        # Sentinel->A
                        if i==0 and j==1 and rank==0:
                            # If this is the first position in the support (support_pos == 0)
                            # and rank > 0 (that is, this is not the first SNP)
                            # and SNPs a, b are adjacent
                            hansel.add_observation('_', snp_a, 0, 1)
                            hansel.add_observation(snp_a, snp_b, 1, 2)
                            crumbs += 1

                        # B->Sentinel
                        elif (j+rank+1) == vcf_handler["N"] and abs(i-j)==1:
                            # Last observation (abs(i-j)==1),
                            # that ends on the final SNP (j+rank+1 == N)
                            hansel.add_observation(snp_a, snp_b, vcf_handler["N"]-1, vcf_handler["N"])
                            hansel.add_observation(snp_b, '_', vcf_handler["N"], vcf_handler["N"]+1)
                            crumbs += 1

                        # A regular observation (A->B)
                        else:
                            hansel.add_observation(snp_a, snp_b, i+rank+1, j+rank+1)
                            crumbs += 1

                            if use_end_sentinels:
                                if j==(support_len-1) and abs(i-j)==1:
                                    # The last SNP on a read, needs a sentinel afterward
                                    hansel.add_observation(snp_b, '_', j+rank+1, j+rank+2)

    bam_queue = Queue()
    progress_queue = Queue()

    # Queue the wokers
    # TODO Evenly divide, but in future, consider the distn
    # TODO Also consider in general block0 has more work to do
    window_l = round((end_pos - start_pos) / float(n_threads))
    for window_i, window_pos in enumerate(range(start_pos, end_pos+1, window_l)):
        bam_queue.put({
            "start": window_pos,
            "end": window_pos + window_l - 1, # add -1 to stop end of window colliding with next window
            "i": window_i,
            "region_end": end_pos,
        })

    processes = []
    for _ in range(n_threads):
        p = Process(target=bam_worker,
                    args=(bam_queue, progress_queue, _))
        processes.append(p)

    # ...and a progress process
    n_reads = Value('i', 0)
    n_observations = Value('i', 0)
    total_covered_snps = Value('i', 0)
    p = Process(target=progress_worker,
                args=(progress_queue, n_threads, n_reads, total_covered_snps, n_observations))
    processes.append(p)

    for p in processes:
        p.start()

    # Add sentinels
    for _ in range(n_threads):
        bam_queue.put(None)

    # Wait for processes to complete work
    for p in processes:
        p.join()


    hansel.n_slices = n_reads.value
    hansel.n_crumbs = n_observations.value
    sys.stderr.write("[NOTE] Loaded %d breadcrumbs from %d bread slices.\n" % (hansel.n_crumbs, hansel.n_slices))

    hansel.L = int(ceil(float(total_covered_snps.value)/n_reads.value))
    sys.stderr.write("[NOTE] Setting Gretel.L to %d\n" % hansel.L)
    return hansel

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


def process_vcf(vcf_path, contig_name, start_pos, end_pos):
    """
    Parse a VCF to extract the genomic positions of called variants.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.

    contig_name : str
        Name of the target contig on which variants were called.

    start_pos : int
        The 1-indexed genomic position from which to begin considering variants.

    end_pos : int
        The 1-indexed genomic position at which to stop considering variants.

    Returns
    -------
    Gretel Metastructure : dict
        A collection of structures used for the execution of Gretel.
        The currently used keys are:
            N : int
                The number of observed SNPs
            snp_fwd : dict{int, int}
                A reverse lookup from the n'th variant, to its genomic position on the contig
            snp_rev : dict{int, int}
                A forward lookup to translate the n'th genomic position to its i'th SNP rank
            region : list{int}
                A masked representation of the target contig, positive values are variant positions
    """

    # Open the VCF
    fp = open(vcf_path, 'rb') # assumes bgzip and tabix
    vcf_records = vcf.Reader(fp)
    n_snps = 0
    snp_reverse = {}
    snp_forward = {}
    region = np.zeros(end_pos + 1, dtype=int)
    i = 0
    for record in vcf_records.fetch(contig_name, 0, end_pos):
        if record.POS < start_pos:
            continue
        if record.POS > end_pos:
            continue

        n_snps += 1
        region[record.POS] = 1
        snp_reverse[i] = record.POS
        snp_forward[record.POS] = i
        i += 1
    fp.close()

    return {
        "N": n_snps,
        "snp_fwd": snp_forward,
        "snp_rev": snp_reverse,
        "region": region,
    }
