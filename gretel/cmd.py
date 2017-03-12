import argparse
import numpy as np
import sys

import gretel
import util

def main():
    """Gretel: A metagenomic haplotyper."""
    parser = argparse.ArgumentParser(description="Gretel: A metagenomic haplotyper.")
    parser.add_argument("bam")
    parser.add_argument("vcf")
    parser.add_argument("contig")
    parser.add_argument("-s", "--start", type=int, default=1, help="1-indexed start base position [default: 1]")
    parser.add_argument("-e", "--end", type=int, default=-1, help="1-indexed end base position [default: contig end]")

    parser.add_argument("-l", "--lorder", type=int, default=0, help="Order of markov chain to predict next nucleotide [default:1]")
    parser.add_argument("-p", "--paths", type=int, default=100, help="Maximum number of paths to generate [default:100]")

    parser.add_argument("--master", default=None, help="Master sequence if available (required to generate out.fasta)")

    parser.add_argument("--quiet", default=False, action='store_true', help="Don't output anything other than a single summary line.")
    parser.add_argument("--sentinels", default=False, action='store_true', help="Add additional sentinels for read ends [default:False][EXPERIMENTAL]")
    parser.add_argument("-o", "--out", default=".", help="Output directory [default .]")
    parser.add_argument("-@", "--threads", type=int, default=1, help="Number of BAM iterators [default 1]")

    ARGS = parser.parse_args()

    VCF_h = gretel.process_vcf(ARGS.vcf, ARGS.contig, ARGS.start, ARGS.end)
    BAM_h = gretel.process_bam(VCF_h, ARGS.bam, ARGS.contig, ARGS.start, ARGS.end, ARGS.lorder, ARGS.sentinels, ARGS.threads)

    #print "[META] #CONTIG", ARGS.contig
    #print "[META] #SNPS", VCF_h["N"]
    #print "[META] #READS", BAM_h["N"]

    PATHS = []
    PATH_PROBS = []
    PATH_PROBS_UW = []
    PATH_FALLS = []

    # Spew out exciting information about the SNPs
    all_marginals = {
        "A": [],
        "C": [],
        "G": [],
        "T": [],
        "N": [],
        "-": [],
        "_": [],
        "total": [],
    }
    if not ARGS.quiet:
        print "i\tpos\tgap\tA\tC\tG\tT\tN\t-\t_\ttot"
        last_rev = 0
        for i in range(0, VCF_h["N"]+1):
            marginal = BAM_h["read_support"].get_counts_at(i)
            snp_rev = 0
            if i > 0:
                snp_rev = VCF_h["snp_rev"][i-1]
            print "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (
                i,
                snp_rev,
                snp_rev - last_rev,
                marginal.get("A", 0),
                marginal.get("C", 0),
                marginal.get("G", 0),
                marginal.get("T", 0),
                marginal.get("N", 0),
                marginal.get("-", 0),
                marginal.get("_", 0),
                marginal.get("total", 0),
            )
            all_marginals["A"].append(marginal.get("A", 0))
            all_marginals["C"].append(marginal.get("C", 0))
            all_marginals["G"].append(marginal.get("G", 0))
            all_marginals["T"].append(marginal.get("T", 0))
            all_marginals["N"].append(marginal.get("N", 0))
            all_marginals["-"].append(marginal.get("-", 0))
            all_marginals["_"].append(marginal.get("_", 0))
            all_marginals["total"].append(
                marginal.get("total", 0)
            )
            last_rev = snp_rev


    # Make some genes
    SPINS = ARGS.paths
    ongoing_mag = 0
    for i in range(0, SPINS):
        init_path, init_prob, init_min = gretel.generate_path(VCF_h["N"], BAM_h["read_support"], BAM_h["read_support_o"])
        if init_path == None:
            break
        current_path = init_path

        MIN_REMOVE = 0.01 # 1%
        if init_min < MIN_REMOVE:
            sys.stderr.write("[RWGT] Ratio %.10f too small, adjusting to %.3f\n" % (init_min, MIN_REMOVE))
            init_min = MIN_REMOVE
        rw_magnitude = gretel.reweight_hansel_from_path(BAM_h["read_support"], init_path, init_min)

        #TODO Horribly inefficient.
        if current_path in PATHS:
            continue
        else:
            ongoing_mag += rw_magnitude
            PATHS.append(current_path)
            PATH_PROBS.append(init_prob["weighted"])
            PATH_PROBS_UW.append(init_prob["unweighted"])
            PATH_FALLS.append(ongoing_mag)
            ongoing_mag = 0
        ongoing_mag += rw_magnitude


    # Make some pretty pictures
    dirn = ARGS.out + "/"
    if ARGS.master:
        master_fa = util.load_fasta(ARGS.master)
        master_seq = master_fa.fetch(master_fa.references[0])
        fasta_out_fh = open(dirn+"out.fasta", "w")
        for i, path in enumerate(PATHS):
            seq = list(master_seq[:])
            for j, mallele in enumerate(path[1:]):
                snp_pos_on_master = VCF_h["snp_rev"][j]
                try:
                    if mallele == "-":
                        # It's a deletion, don't print a SNP
                        seq[snp_pos_on_master-1] = ""
                    else:
                        seq[snp_pos_on_master-1] = mallele
                except IndexError:
                    print path, len(seq), snp_pos_on_master-1
                    sys.exit(1)
            fasta_out_fh.write(">%d__%.2f\n" % (i, PATH_PROBS[i]))
            fasta_out_fh.write("%s\n" % "".join(seq[ARGS.start-1 : ARGS.end]))
        fasta_out_fh.close()
    else:
        fasta_out_fh = open(dirn+"out.fasta", "w")
        for i, path in enumerate(PATHS):
            fasta_out_fh.write(">%d__%.2f\n" % (i, PATH_PROBS[i]))
            fasta_out_fh.write("%s\n" % "".join("".join(path[1:])))
        fasta_out_fh.close()

    #TODO datetime, n_obs, n_slices, avg_obs_len, L, n_paths, n_avg_loglik
    crumb_file = open(dirn+"gretel.crumbs", "w")
    crumb_file.write("# %d\t%d\t%d\t%.2f\t%d\t%.2f\t%.2f\t%.2f\n" % (
        VCF_h["N"],
        BAM_h["read_support"].n_crumbs,
        BAM_h["read_support"].n_slices,
        BAM_h["meta"]["L"],
        BAM_h["read_support"].L,
        np.mean(PATH_PROBS),
        np.mean(PATH_PROBS_UW),
        np.mean(PATH_FALLS),
    ))
    for p in range(len(PATHS)):
        crumb_file.write("%d\t%.2f\t%.2f\t%.2f\n" % (
                p,
                PATH_PROBS[p],
                PATH_PROBS_UW[p],
                PATH_FALLS[p]
        ))
