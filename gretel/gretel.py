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
    size = 0
    for i in range(0, n_snps+1):
        for j in range(0, i+1+1):
            # Reduce read supports
            if i >= n_snps:
                size += supports_mat.reweight_observation(path[i], path[j], i, i+1, ratio)
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
                size += supports_mat.reweight_observation(path[t_i], path[t_j], t_i, t_j, ratio)
    return size


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

    meta = util.load_from_bam(hansel, bam, contig_name, vcf_handler)

    return {
        "read_support": hansel,
        "read_support_o": hansel.copy(),
        "meta": meta,
    }

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

