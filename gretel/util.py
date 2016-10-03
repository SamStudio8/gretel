import pysam
import numpy as np
import sys

#TODO SENTINEL SYMBOLS BEFORE AND AFTER A READ
#TODO What happens if we traverse backwards...?
#TODO Single SNP reads could use a pairwise observation with themselves? (A, A, i, i)

def load_from_bam(hansel, bam, target_contig, start_pos, end_pos, vcf_handler, use_end_sentinels=False):
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
    support_seq_sum = 0
    support_seq_num = 0
    for read in bam.fetch(target_contig):
        START_POS_OFFSET = 0

        hansel.n_slices += 1

        if read.is_duplicate or read.is_secondary:
            continue

        LEFTMOST_1pos = read.reference_start + 1 # Convert 0-based reference_start to 1-based position (to match region array and 1-based VCF)
        if LEFTMOST_1pos < start_pos:
            # Read starts before the start_pos
            if read.reference_start + 1 + read.query_alignment_length < start_pos:
                # Read ends before the start_pos
                continue
            LEFTMOST_1pos = start_pos
            START_POS_OFFSET = (start_pos - (read.reference_start + 1))
        elif LEFTMOST_1pos > end_pos:
            # The BAM is sorted and this read begins after the end of the region of interest...
            break

        #RIGHTMOST_1pos = read.reference_start + 1 + read.query_alignment_length
        RIGHTMOST_1pos = read.reference_end #ofc this is 1-indexed instead of 0
        if RIGHTMOST_1pos > end_pos:
            # Read ends after the end_pos
            RIGHTMOST_1pos = end_pos

        # Check there is any support
        support_len = np.sum(vcf_handler["region"][LEFTMOST_1pos : RIGHTMOST_1pos + 1])

        # Ignore reads without evidence
        if support_len == 0:
            continue

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

        support_seq_num += 1
        support_seq_sum += len(support_seq.replace("N", "").replace("-", ""))

        print read.reference_start + 1

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
                if snp_a in hansel.unsymbols:
                    continue

                # Sentinel->A
                if i==0 and j==1 and rank==0:
                    # If this is the first position in the support (support_pos == 0)
                    # and rank > 0 (that is, this is not the first SNP)
                    # and SNPs a, b are adjacent
                    hansel.add_observation('_', snp_a, 0, 1)
                    hansel.add_observation(snp_a, snp_b, 1, 2)

                # B->Sentinel
                elif (j+rank+1) == vcf_handler["N"] and abs(i-j)==1:
                    # Last observation (abs(i-j)==1),
                    # that ends on the final SNP (j+rank+1 == N)
                    hansel.add_observation(snp_a, snp_b, vcf_handler["N"]-1, vcf_handler["N"])
                    hansel.add_observation(snp_b, '_', vcf_handler["N"], vcf_handler["N"]+1)

                # A regular observation (A->B)
                else:
                    hansel.add_observation(snp_a, snp_b, i+rank+1, j+rank+1)

                    if use_end_sentinels:
                        if j==(support_len-1) and abs(i-j)==1:
                            # The last SNP on a read, needs a sentinel afterward
                            hansel.add_observation(snp_b, "_", j+rank+1, j+rank+2)


    meta["support_seq_avg"] = support_seq_sum/float(support_seq_num)
    if hansel.L == 0:
        from math import ceil
        hansel.L = int(ceil(meta["support_seq_avg"])) #TODO
        sys.stderr.write("[NOTE] Setting Gretel.L to %d\n" % hansel.L)
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

