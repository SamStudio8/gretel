import argparse
import matplotlib.pyplot as plt
import numpy as np

from hansel import Hansel
import gretel

def main():
    parser = argparse.ArgumentParser(description="Gretel: A metagenomic haplotyper.")
    parser.add_argument("bam")
    parser.add_argument("vcf")
    parser.add_argument("contig")
    parser.add_argument("-s", "--start", type=int, default=1, help="1-indexed start base position [default: 1]")
    parser.add_argument("-e", "--end", type=int, default=-1, help="1-indexed end base position [default: contig end]")

    parser.add_argument("--genes", default=None, help="Input genes for verification")
    parser.add_argument("--hit", default=None, help="Hit table for verification")

    ARGS = parser.parse_args()

    # Process hit table and FASTA reference (if provided)
    HITS = []
    REFS = []
    REF_NAMES = []
    if ARGS.hit and ARGS.genes:
        HITS = gretel.process_hits(ARGS.hit)
        REFS, REF_NAMES = gretel.process_refs(ARGS.genes)


    VCF_h = gretel.process_vcf(ARGS.vcf, ARGS.contig, ARGS.start, ARGS.end)
    BAM_h = gretel.process_bam(VCF_h, ARGS.bam, ARGS.contig)

    print "[META] #CONTIG", ARGS.contig
    print "[META] #SNPS", VCF_h["N"]
    #print "[META] #READS", BAM_h["N"]

    PATHS = []
    PATH_PROBS = []
    PATH_PROBS_UW = []

    # Spew out exciting information about the SNPs
    all_marginals = {
        "A": [],
        "C": [],
        "G": [],
        "T": [],
        "N": [],
        "total": [],
    }
    print "i\tA\tC\tG\tT\tN"
    for i in range(0, VCF_h["N"]+1):
        marginal = BAM_h["read_support"].get_counts_at(i)
        print "%d\t%d\t%d\t%d\t%d\t%d" % (
            i,
            marginal.get("A", 0),
            marginal.get("C", 0),
            marginal.get("G", 0),
            marginal.get("T", 0),
            marginal.get("N", 0),
        )
        all_marginals["A"].append(marginal.get("A", 0))
        all_marginals["C"].append(marginal.get("C", 0))
        all_marginals["G"].append(marginal.get("G", 0))
        all_marginals["T"].append(marginal.get("T", 0))
        all_marginals["N"].append(marginal.get("N", 0))
        all_marginals["total"].append(
            marginal.get("total", 0)
        )


    # Make some genes
    SPINS = 100
    for i in range(0, SPINS):
        init_path, init_prob, init_min = gretel.establish_path(VCF_h["N"], BAM_h["read_support"], BAM_h["read_support_o"])
        if init_path == None:
            break
        current_path = init_path
        BAM_h["read_support"] = gretel.add_ignore_support3(BAM_h["read_support"], VCF_h["N"], init_path, init_min)
        PATHS.append(current_path)
        PATH_PROBS.append(init_prob["weighted"])
        PATH_PROBS_UW.append(init_prob["unweighted"])


    # Make some pretty pictures
    import matplotlib.pyplot as plt

    fasta_out_fh = open("out.fasta", "w")
    master_fh = open("master.fa")
    master_fh.readline()
    master_seq = "".join([l.strip() for l in master_fh.readlines()])
    master_fh.close()
    for i, path in enumerate(PATHS):
        seq = list(master_seq)
        for j, mallele in enumerate(path[1:]):
            snp_pos_on_master = VCF_h["snp_rev"][j]
            seq[snp_pos_on_master-1] = mallele
        fasta_out_fh.write(">%d__%.2f\n" % (i, PATH_PROBS[i]))
        fasta_out_fh.write("%s\n" % "".join(seq))
    fasta_out_fh.close()

    #TODO None for ARGS
    con_mat, snp_mat, seen_mat, master_mat, full_confusion = gretel.confusion_matrix(PATHS, VCF_h, HITS, REFS, REF_NAMES, VCF_h["N"])

    print "#\tname\tsites\trate\tbestit0\trefd\tlogl\tmap"
    RECOVERIES = []
    for i in range(0, len(con_mat)):
        recovered = 0
        at = None
        for j in range(0, len(con_mat[i])):
            if con_mat[i][j] > recovered:
                recovered = con_mat[i][j]
                at = j
        if at != None:
            RECOVERIES.append(at)
            print "%d\t%s\t%d\t%.2f\t%d\t%d\t%.2f\t%s" % (
                    i,
                    REF_NAMES[i][1:31],
                    np.mean(seen_mat[i]),
                    recovered,
                    at,
                    np.sum(master_mat[i][at]),
                    PATH_PROBS_UW[at],
                    "".join([str(int(n)) for n in full_confusion[i][at]]),
            )

    fig, ax = plt.subplots(3,1,sharex=True)
    x_ax = range(0, len(PATHS))
    ax[0].set_title("Gene Identity Recovery by Iteration")
    for i in range(0, len(REFS)):
        ax[0].plot(x_ax, con_mat[i], linewidth=2.0, alpha=0.75)
    ax[0].set_ylabel("Identity (%)")
    ax[0].set_ylim(0, 100)

    PATH_PROBS_RATIO = []
    for i, p in enumerate(PATH_PROBS):
        try:
            PATH_PROBS_RATIO.append( PATH_PROBS[i] / PATH_PROBS_UW[i] )
        except:
            PATH_PROBS_RATIO.append(0)
    # Add likelihood
    ax[1].plot(x_ax, PATH_PROBS, color="red", linewidth=3.0)
    ax[1].plot(x_ax, PATH_PROBS_UW, color="green", linewidth=3.0)
    ax[1].set_ylabel("Log(P)")
    ax[1].set_title("Path Likelihood by Iteration")

    ax[2].plot(x_ax, PATH_PROBS_RATIO, color="blue", linewidth=3.0)
    ax[2].set_xlabel("Iteration (#)")

    # Add recoveries
    for r in RECOVERIES:
        ax[0].axvline(r, color='k', linestyle='--')
        ax[1].axvline(r, color='k', linestyle='--')
    plt.show()

    """
    for i in range(0, len(REFS)):
        plt.pcolor(snp_mat[i], cmap=plt.cm.Blues)
        plt.show()
    for i in range(0, len(PATHS)):
        plt.pcolor(snp_mat[i], cmap=plt.cm.Blues)
        print snp_mat[i]
        plt.show()

    running_bottom = None
    plt.bar(range(0, VCF_h["N"]+1), np.array(all_marginals["A"])/np.array(all_marginals["total"]), color="blue")
    running_bottom = np.array(all_marginals["A"])/np.array(all_marginals["total"])

    plt.bar(range(0, VCF_h["N"]+1), np.array(all_marginals["C"])/np.array(all_marginals["total"]), bottom=running_bottom, color="green")
    running_bottom += np.array(all_marginals["C"])/np.array(all_marginals["total"])

    plt.bar(range(0, VCF_h["N"]+1), np.array(all_marginals["G"])/np.array(all_marginals["total"]), bottom=running_bottom, color="red")
    running_bottom += np.array(all_marginals["G"])/np.array(all_marginals["total"])

    plt.bar(range(0, VCF_h["N"]+1), np.array(all_marginals["T"])/np.array(all_marginals["total"]), bottom=running_bottom, color="yellow")
    running_bottom += np.array(all_marginals["T"])/np.array(all_marginals["total"])

    plt.bar(range(0, VCF_h["N"]+1), np.array(all_marginals["N"])/np.array(all_marginals["total"]), bottom=running_bottom, color="black")
    plt.show()
    """
