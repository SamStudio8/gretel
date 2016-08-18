import pysam
import numpy as np
import sys

#TODO SENTINEL SYMBOLS BEFORE AND AFTER A READ
#TODO What happens if we traverse backwards...?
#TODO Single SNP reads could use a pairwise observation with themselves? (A, A, i, i)

def load_from_bam(hansel, bam, target_contig, vcf_handler, use_end_sentinels=False):
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
    support_seq_lens = []
    for read in bam.fetch(target_contig):
        hansel.n_slices += 1

        if read.is_duplicate or read.is_secondary:
            continue

        # Check there is any support
        LEFTMOST_1pos = read.reference_start + 1 # Convert 0-based reference_start to 1-based position (to match region array and 1-based VCF)
        support_len = sum(vcf_handler["region"][LEFTMOST_1pos: LEFTMOST_1pos+read.query_alignment_length])

        # Ignore reads without evidence
        if support_len == 0:
            continue

        rank = sum(vcf_handler["region"][0 : LEFTMOST_1pos])

        support_seq = ""
        for i in range(0, support_len):
            offset = 0
            snp_rev = vcf_handler["snp_rev"][rank + i]
            snp_pos_on_read = snp_rev - LEFTMOST_1pos

            bases_observed = 0          # Bases observed via CIGAR so far
            last_offsetting_op = 0      # Did we under or overshoot the target?

            for cigar in read.cigartuples:
                # Current operation type and number of bases
                op, count = cigar
                bases_observed += count

                if op == 0:
                    # Match
                    pass
                elif op == 1:
                    # Insert
                    #offset += count  # Bases appearing 'later' than expected
                    #last_offsetting_op = op
                    pass
                elif op == 2:
                    # Deletion
                    #offset -= count  # Bases appearing 'earlier' than expected
                    #last_offsetting_op = op
                    pass
                elif op == 4:
                    # Soft Clip
                    pass
                else:
                    raise Exception("Unsupported CIGAR Opcode (%d) Encountered on '%s'" % (op, read.qname))

                if bases_observed >= snp_pos_on_read:  # TODO >= ?
                    # Abort early if we find the target SNP site
                    break

            # We should have overshot
            if bases_observed >= snp_pos_on_read:
                pass
            else:
                raise Exception("Failed to reach a SNP site (%d) on '%s'" % (snp_rev, read.qname))

            support_seq += read.query_alignment_sequence[snp_pos_on_read + offset]

        support_seq_lens.append(len(support_seq.replace("N", "").replace("-", "")))

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


    meta["support_seq_lens"] = support_seq_lens
    if hansel.L == 0:
        from math import ceil
        hansel.L = int(ceil(np.mean(support_seq_lens))) #TODO
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

