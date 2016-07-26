import sys
from math import log,log10,exp
import random

import vcf
import pysam
import numpy as np

from hansel import Hansel
import util

#TODO Look at some SNPs outside the gene...
#TODO Reimplement SNP ladder generator with new matrix format
#TODO Should the denom of the conditional use the unique variants at i-l or i?
#TODO Util to parse known input and return SNP seq

## PROBABILITY ################################################################
def add_ignore_support3(supports_mat, n_snps, path, ratio):
    for i in range(0, n_snps+1):
        for j in range(0, i+1+1):
            # Reduce read supports
            if i >= n_snps:
                supports_mat.reweight_observation(path[i], path[j], i, i+1, ratio)
                break #???
            else:
                if j < i:
                    # This isn't just a case of j < i, but means that we are looking
                    # at the two SNPs the wrong way around, we must switch them before
                    # handing them over to reweight_observation
                    t_i = j
                    t_j = i
                else:
                    t_i = i
                    t_j = j
                supports_mat.reweight_observation(path[t_i], path[t_j], t_i, t_j, ratio)


## INPUT OUTPUT ###############################################################
def process_hits(hit_tab):
    HITS = []
    with open(hit_tab) as hit_fh:
        for line in hit_fh:
            if line[0] == "#":
                continue

            line = line.strip()
            if len(line) == 0:
                continue

            fields = line.split("\t")
            HITS.append({
                "subject": fields[1],
                "length": int(fields[3]),
                "ref_s": int(fields[6]),
                "ref_e": int(fields[7]),
                "sub_s": int(fields[8]),
                "sub_e": int(fields[9]),
            })
    return HITS


def process_refs(ref_fa):
    return util.load_fasta(ref_fa)

def process_vcf(vcf_path, contig_name, start_pos, end_pos):
    # Open the VCF
    vcf_records = vcf.Reader(open(vcf_path))
    n_snps = 0
    snp_reverse = {}
    snp_forward = {}
    region = np.zeros(end_pos - start_pos + 2, dtype=int)
    for i, record in enumerate(vcf_records.fetch(contig_name, start_pos-1, end_pos)):
        n_snps += 1
        region[record.POS] = 1
        snp_reverse[i] = record.POS
        snp_forward[record.POS] = i

    return {
        "N": n_snps,
        "snp_fwd": snp_forward,
        "snp_rev": snp_reverse,
        "region": region,
    }

def process_bam(vcf_handler, bam_path, contig_name, L):
    bam = pysam.AlignmentFile(bam_path)

    #NOTE(samstudio8)
    # Could we optimise for lower triangle by collapsing one of the dimensions
    # such that Z[m][n][i][j] == Z[m][n][i + ((j-1)*(j))/2]
    read_support_mat = np.zeros( (6, 6, vcf_handler["N"]+2, vcf_handler["N"]+2) )
    hansel = Hansel(read_support_mat, ['A', 'C', 'G', 'T', 'N', "_"], ['N', "_"], L=L)

    util.load_from_bam(hansel, bam, contig_name, vcf_handler)

    return {
        "read_support": hansel,
        "read_support_o": hansel.copy(),
    }

def confusion_matrix(PATHS, VCF_h, HITS, REFS, REF_NAMES, N, master_path=None):
    #FIXME Kinda gross
    CONFUSION = np.zeros( (len(REFS), len(PATHS)) )
    SEEN = np.zeros( (len(REFS), len(PATHS)) )
    snp_matrix = np.zeros( (len(REFS), len(PATHS), N) )
    mat_matrix = np.zeros( (len(REFS), len(PATHS), N) )
    #snp_matrix = np.zeros( (len(PATHS), len(REFS), N+1) )
    full_confusion = np.zeros(( len(REFS), len(PATHS), N ))

    if master_path:
        master_fa = util.load_fasta(master_path)
        master_seq = master_fa.fetch(master_fa.references[0])

    if not (HITS and REFS):
        sys.stderr.write("Cannot call confusion_matrix without references to check against :<\n")

    for path_id, path_variants in enumerate(PATHS):
        for ref_gene_id, ref_gene_name in enumerate(REF_NAMES):

            # Get the current input gene sequence
            ref_seq = REFS.fetch(ref_gene_name)

            for hit_record in [h for h in HITS if h["subject"] == ref_gene_name]:
                # Test a path, against a hit for the current reference
                # Typically there is just one hit for each reference, but this is not strictly true.
                sys.stderr.write("[TEST] PATH%d, %s with hit %s\n" % (path_id, ref_gene_name, str(hit_record)))

                # Iterate over the current path and compare the malleles
                # (Ignore the first variant of the path, it's the sentinel)
                for variant_id, variant in enumerate(path_variants[1:]):
                    full_confusion[ref_gene_id][path_id][variant_id] = 9

                    # Translate the variant_id to 1-based position of the variant on the MASTER (pseudo-ref)
                    snp_pos_on_master = VCF_h["snp_rev"][variant_id]

                    # Does the hit cover the reference at the variant?
                    if snp_pos_on_master < hit_record["ref_s"] or snp_pos_on_master > hit_record["ref_e"]:
                        continue

                    #TODO should probably consider ARGS["REGION_START"]
                    # Translate 1-based MASTER position to 0-based (ref_seq string) based on hit_record:
                    # ...convert MASTER position (1-indexed) to 0-indexed (ref string is python string)
                    position = snp_pos_on_master - 1

                    # ...add offset of subject start (-1 to convert from blast hit6 1-pos to 0-pos)
                    #   The subject has some leader sequence that does not hit the
                    #   query (reference), we ignore it by adding the sub_s offset.
                    #
                    #     >>>>>>>>>>>>>>|================================ REF(QRY)
                    #     |============================================== HIT(SUB)
                    #                   |
                    #                   *sub_s (1-indexed start of hit on subject)
                    position += hit_record["sub_s"] - 1

                    # ...remove offset of reference start (-1 to convert from blast hit6 1-pos to 0-pos)
                    #   The subject does not cover the entirity of the query (reference),
                    #   positions of the subject therefore refer to earlier parts of
                    #   the query, we adjust this by shifting positions backward...
                    #
                    #              *ref_s (1-indexed start of hit on query)
                    #              |
                    #     |============================================== REF(QRY)
                    #     <<<<<<<<<|===================================== HIT(SUB)
                    position -= hit_record["ref_s"] - 1

                    # Does the variant at the reference site match the variant of the current path?
                    if ref_seq[position] == variant:
                        if master_path:
                            if variant != master_seq[snp_pos_on_master-1]:
                                mat_matrix[ref_gene_id][path_id][variant_id] = 1
                        CONFUSION[ref_gene_id][path_id] += 1
                        full_confusion[ref_gene_id][path_id][variant_id] = 1
                        snp_matrix[ref_gene_id][path_id][variant_id] = 100
                        #snp_matrix[pi][ri][i] += 100
                    else:
                        full_confusion[ref_gene_id][path_id][variant_id] = 0
                    #snp_matrix[pi][ri][i] += 50
                    SEEN[ref_gene_id][path_id] += 1

    return (CONFUSION / SEEN) * 100, snp_matrix, SEEN, mat_matrix, full_confusion

## PATH GENERATION ############################################################

def append_path(path, next_m, next_v):
    # Probably a bit gross as it has side effects on path...
    if next_m is not None:
        path.append(next_m)
    else:
        raise Exception("Cowardly refusing to append None as a nucleotide. Cheerio.")

def establish_path(n_snps, read_support_mat, original_read_support):
    # Cross the metahaplome in a greedily, naive fashion to establish a base path
    # This seeds the rest of the path generation (we might want to just select
    #   a random path here in future)

    running_prob = 0.0
    running_prob_uw = 0.0
    current_path = ['_'] # start with the dummy
    marginals = []
    L = 5

    # Find path
    sys.stderr.write("*** ESTABLISH ***\n")
    for snp in range(1, n_snps+1):
        sys.stderr.write("\t*** ***\n")
        sys.stderr.write("\t[SNP_] SNP %d\n" % snp)

        # Get marginal and calculate branch probabilities for each available
        # mallele, given the current path seen so far
        # Select the next branch and append it to the path
        curr_branches = read_support_mat.get_edge_weights_at(snp, current_path)
        sys.stderr.write("\t[TREE] %s\n" % curr_branches)
        # Return the symbol and probability of the next base to add to the
        # current path based on the best marginal
        next_v = 0.0
        next_m = None

        for symbol in curr_branches:
            if symbol == "total":
                continue
            if next_m is None:
                next_v = curr_branches[symbol]
                next_m = symbol
            elif curr_branches[symbol] > next_v:
                next_v = curr_branches[symbol]
                next_m = symbol

        if next_m == None:
            sys.stderr.write("[FAIL] Unable to select next branch from %d to %d\n" % (snp-1, snp))
            return None, None, None
            current_marg = read_support_mat.get_counts_at(snp)
            next_m = random.choice(['A', 'C', 'G', 'T'])
            false_marg_ratio = 1 / (1+current_marg["total"])

            marginals.append(false_marg_ratio)
            running_prob += log10(false_marg_ratio)
            current_marg = original_read_support.get_counts_at(snp)
            false_marg_ratio = 1 / (1+current_marg["total"])
            running_prob_uw += log10(false_marg_ratio)
#TODO testing
        else:
            selected_edge_weight = read_support_mat.get_marginal_of_at(next_m, snp)
            marginals.append(selected_edge_weight) #TODO This isn't a log, is it accurate enough later?

            running_prob += log10(selected_edge_weight)
            running_prob_uw += log10(original_read_support.get_marginal_of_at(next_m, snp))
        append_path(current_path, next_m, next_v)

    return current_path, {"unweighted": running_prob_uw, "weighted": running_prob}, min(marginals)

