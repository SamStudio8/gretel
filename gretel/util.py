import pysam
import numpy as np
import sys

#TODO SENTINEL SYMBOLS BEFORE AND AFTER A READ
#TODO What happens if we traverse backwards...?

def load_from_bam(hansel, bam, target_contig, vcf_handler):
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

                    if j==(support_len-1) and abs(i-j)==1:
                        # The last SNP on a read, needs a sentinel afterward
                        #hansel.add_observation(snp_b, "_", j+rank+1, j+rank+2)
                        pass


    if hansel.L == 1:
        from math import ceil
        hansel.L = int(ceil(np.median(support_seq_lens))) #TODO
        sys.stderr.write("[NOTE] Setting Gretel.L to %d\n" % hansel.L)

def load_fasta(fa_path):
    return pysam.FastaFile(fa_path)

