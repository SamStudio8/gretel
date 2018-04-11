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
    parser.add_argument("-e", "--end", type=int, default=-1, help="1-indexed end base position [default: reference length]")

    parser.add_argument("-l", "--lorder", type=int, default=0, help="Order of markov chain to predict next nucleotide [default calculated from read data]")
    parser.add_argument("-p", "--paths", type=int, default=100, help="Maximum number of paths to generate [default:100]")

    parser.add_argument("--master", default=None, help="Master sequence (will be used to fill in homogeneous gaps in haplotypes, otherwise Ns)") #TODO Use something other than N? Should probably be a valid IUPAC

    parser.add_argument("--quiet", default=False, action='store_true', help="Don't output anything other than a single summary line.")
    parser.add_argument("--sentinels", default=False, action='store_true', help="Add additional sentinels for read ends [default:False][EXPERIMENTAL]")
    parser.add_argument("-o", "--out", default=".", help="Output directory [default .]")
    parser.add_argument("-@", "--threads", type=int, default=1, help="Number of BAM iterators [default 1]")

    ARGS = parser.parse_args()

    if ARGS.end == -1:
        ARGS.end = util.get_ref_len_from_bam(ARGS.bam, ARGS.contig)
        sys.stderr.write("[NOTE] Setting end_pos to %d" % ARGS.end)

    VCF_h = gretel.process_vcf(ARGS.vcf, ARGS.contig, ARGS.start, ARGS.end)
    BAM_h = gretel.process_bam(VCF_h, ARGS.bam, ARGS.contig, ARGS.start, ARGS.end, ARGS.lorder, ARGS.sentinels, ARGS.threads)

    # Check if there is a gap in the matrix
    # TODO(samstudio8) Ideally we would do this IN the bam worker threads such
    #                  that a thread could raise an error and halt the whole
    #                  process if we already know we'll yield a matrix we cannot
    #                  actually work with, but this will do for now...
    for i in range(0, VCF_h["N"]+1):
        marginal = BAM_h["read_support"].get_counts_at(i)

        if i > 0:
            snp_rev = VCF_h["snp_rev"][i-1]
        else:
            snp_rev = 0
        if marginal.get("total", 0) == 0:
            sys.stderr.write('''[FAIL] Unable to recover pairwise evidence concerning SNP #%d at position %d
       Gretel needs every SNP to appear on a read with at least one other SNP, at least once.
       There is no read in your data set that bridges SNP #%d with any of its neighbours.

       * If you are trying to run Gretel along an entire contig or genome, please note that
       this is not the recommended usage for Gretel, as it was intended to uncover the
       variation in a metahaplome: the set of haplotypes for a specific gene.
           See our pre-print https://doi.org/10.1101/223404 for more information

       Consider running a prediction tool such as `prokka` on your assembly or reference
       and using the CDS regions in the GFF for corresponding genes of interest to
       uncover haplotypes with Gretel instead.

       * If you are already doing this, consider calling for SNPs more aggressively.
       We use `snpper` (https://github.com/SamStudio8/gretel-test/blob/master/snpper.py)
       to determine any site in a BAM that has at least one read in disagreement with
       the reference as a SNP. Although this introduces noise from alignment and sequence
       error, Gretel is fairly robust. Importantly, this naive calling method will
       likely close gaps between SNPs and permit recovery.

       * Finally, consider that the gaps are indicative that your reads do not support
       one or more parts of your assembly or reference. You could try and find or construct
       a more suitable reference, or reduce the size of the recovery window.

       Sorry :(\n''' % (i, snp_rev, i))
            sys.exit(1)


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
    fasta_out_fh = open(dirn+"out.fasta", "w")
    if ARGS.master:
        master_fa = util.load_fasta(ARGS.master)
        master_seq = master_fa.fetch(master_fa.references[0])
    else:
        master_seq = ["N"] * ARGS.end

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
