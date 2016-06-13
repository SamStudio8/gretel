import vcf
import pysam
import numpy as np
import sys
from math import log,log10,exp
import random

from hansel import Hansel

#TODO Look at some SNPs outside the gene...
#TODO Reimplement SNP ladder generator with new matrix format
#TODO Should the denom of the conditional use the unique variants at i-l or i?
#TODO Util to parse known input and return SNP seq

## PROBABILITY ################################################################
def add_ignore_support3(supports_mat, n_snps, path, ratio):
    print "********** REWEIGHTING **********"
    print path

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
    return supports_mat


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
    REFERENCES = {}
    with open(ref_fa) as fasta_fh:
        current_seq = []
        current_name = None
        names = []
        for line in fasta_fh:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    REFERENCES[current_name] = "".join(current_seq)

                current_name = line
                names.append(current_name)
                current_seq = []
            else:
                current_seq.append(line)

    REFERENCES[current_name] = "".join(current_seq)
    return REFERENCES, names

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

def process_bam(vcf_handler, bam_path, contig_name):
    bam = pysam.AlignmentFile(bam_path)

    #NOTE(samstudio8)
    # Could we optimise for lower triangle by collapsing one of the dimensions
    # such that Z[m][n][i][j] == Z[m][n][i + ((j-1)*(j))/2]
    read_support_mat = np.zeros( (5, 5, vcf_handler["N"]+2, vcf_handler["N"]+2) )
    hansel = Hansel(read_support_mat)

    hansel.load_from_bam(bam, contig_name, vcf_handler)

    return {
        "read_support": hansel,
        "read_support_o": hansel.copy(),
    }

def confusion_matrix(PATHS, VCF_h, HITS, REFS, REF_NAMES, N):
    #FIXME Kinda gross
    CONFUSION = np.zeros( (len(REFS), len(PATHS)) )
    SEEN = np.zeros( (len(REFS), len(PATHS)) )
    snp_matrix = np.zeros( (len(REFS), len(PATHS), N) )
    mat_matrix = np.zeros( (len(REFS), len(PATHS), N) )
    #snp_matrix = np.zeros( (len(PATHS), len(REFS), N+1) )
    full_confusion = np.zeros(( len(REFS), len(PATHS), N ))

    master_fh = open("master.fa")
    master_fh.readline()
    master_seq = "".join([l.strip() for l in master_fh.readlines()])
    master_fh.close()

    # For each new path
    for pi, PATH in enumerate(PATHS):

        # ...and that we have a hit table with references
        if HITS and REFS:

            # ... go over each reference available
            for ri, reference in enumerate(REF_NAMES):
                # ...and get each hit to that reference
                for hit in [h for h in HITS if h["subject"] in reference]:
                    print "[TEST] PATH%d, %s with hit %s" % (pi, reference, str(hit))
                    # Get the corresponding part of the reference
                    #ref_ = REFS[reference][hit["ref_s"]-1:hit["ref_e"]]
                    ref_ = REFS[reference]

                    # Get the corresponding part of the new path
                    # Check the corresponding SNPs of the new path against ref_
                    for i, mallele in enumerate(PATH[1:]):
                        full_confusion[ri][pi][i] = 9

                        # Does the hit cover the reference at the SNP?
                        snp_pos_on_master = VCF_h["snp_rev"][i]
                        if snp_pos_on_master < hit["ref_s"] or snp_pos_on_master > hit["ref_e"]:
                            continue


                        #TODO should probably consider ARGS["REGION_START"]
                        position = (hit["sub_s"]-1) + snp_pos_on_master-1 - (hit["ref_s"]-1)

                        #print i, mallele, VCF_h["snp_rev"][i], ARGS["REGION_START"], hit["sub_s"]
                        #position = VCF_h["snp_rev"][i]-ARGS["REGION_START"]+hit["sub_s"]-1

                        #if position <= len(ref_) and position >= hit["sub_s"]-1 and position <= hit["sub_e"]:
                        #if position < len(ref_):
                        print i, snp_pos_on_master, position, ref_[position], mallele, ref_[position] == mallele

                        if ref_[position] == mallele:
                            if mallele != master_seq[snp_pos_on_master-1]:
                                mat_matrix[ri][pi][i] = 1
                            CONFUSION[ri][pi] += 1
                            full_confusion[ri][pi][i] = 1
                            snp_matrix[ri][pi][i] = 100
                            #snp_matrix[pi][ri][i] += 100
                        else:
                            full_confusion[ri][pi][i] = 0
                        #snp_matrix[pi][ri][i] += 50
                        SEEN[ri][pi] += 1
            print CONFUSION
            print "**"
            print SEEN
            print mat_matrix
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
    current_path = ['N'] # start with the dummy
    marginals = []
    L = 5

    # Find path
    print "*** ESTABLISH ***"
    for snp in range(1, n_snps+1):
        print "\t*** ***"
        print "\t[SNP_] SNP %d" % snp

        # Get marginal and calculate branch probabilities for each available
        # mallele, given the current path seen so far
        # Select the next branch and append it to the path
        next_m, next_v = read_support_mat.select_next_edge_at(snp, current_path, L=L)
        if next_m == None:
            print "[FAIL] Unable to select next branch from %d to %d" % (snp-1, snp)
            return None, None, None
            current_marg = read_support_mat.get_marginal_at(snp)
            next_m = random.choice(['A', 'C', 'G', 'T'])
            false_marg_ratio = 1 / (1+current_marg["total"])

            marginals.append(false_marg_ratio)
            running_prob += log10(false_marg_ratio)
            current_marg = original_read_support.get_marginal_at(snp)
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

