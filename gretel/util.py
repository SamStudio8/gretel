import pysam

def load_from_bam(hansel, bam, target_contig, vcf_handler):
    for read in bam.fetch(target_contig):
        hansel.n_slices += 1

        if read.is_duplicate or read.is_secondary:
            continue

        # Check there is any support
        support_len = sum(vcf_handler["region"][read.pos+1: read.pos+read.qlen+1])
        if support_len == 0:
            continue

        rank = sum(vcf_handler["region"][0 : read.pos+1])

        support_seq = ""
        for i in range(0, support_len):
            offset = 0
            snp_rev = vcf_handler["snp_rev"][rank + i]
            snp_pos_on_read = snp_rev - read.pos - 1

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

        #TODO ew.
        for i in range(0, support_len):
            snp_a = support_seq[i]
            for j in range(i+rank+1, rank+support_len):  #TODO Start from i+rank ?
                snp_b = support_seq[j-rank]

                if i == 0 and rank > 0 and abs(i+rank-j)==1:
                    hansel.add_observation(snp_a, snp_b, rank+1, rank+2)
                elif j-1 == 0:
                    hansel.add_observation('N', snp_a, 0, 1)
                    hansel.add_observation(snp_a, snp_b, i+rank+1, j+1)
                else:
                    hansel.add_observation(snp_a, snp_b, i+rank+1, j+1)

                if (j+1) == vcf_handler["N"] and abs(i+rank-j) == 1:
                    hansel.add_observation(snp_b, 'N', vcf_handler["N"]-1+1, vcf_handler["N"]+1)

def load_fasta(fa_path):
    return pysam.FastaFile(fa_path)
