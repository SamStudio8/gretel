import sys
import matplotlib.pyplot as plt
hits = open(sys.argv[1])
crumbs = open(sys.argv[2])

PROBS = []
PROBS_UW = []
WEIGHTS = []
for line in crumbs:
    fields = [float(x) for x in line.strip().split("\t")]
    PROBS.append(fields[1])
    PROBS_UW.append(fields[2])
    WEIGHTS.append(fields[3])

mismatches = {}
best_by_id = {}

for line in hits:
    fields = line.strip().split("\t")
    query = fields[1]
    gene = fields[0]
    pid = int(query.split("__")[0])

    if gene not in mismatches:
        mismatches[gene] = [None] * len(PROBS)
    mismatches[gene][pid] = int(fields[4])

    metric = float(fields[11])
    if gene not in best_by_id:
        best_by_id[gene] = {
            "p": query,
            "id": float(fields[2]),
            "line": line.strip(),
            "p_": pid,
            "metric": metric
        }
    else:
        if metric > best_by_id[gene]["metric"]:
            best_by_id[gene] = {
                "p": query,
                "id": float(fields[2]),
                "line": line.strip(),
                "p_": pid,
                "metric": metric
            }


fig, ax = plt.subplots(3,1,sharex=True)
x_ax = range(0, len(mismatches[mismatches.keys()[0]]) )
ax[0].set_title("Gene Mismatches by Iteration")
for g in mismatches:
    ax[0].plot(x_ax, mismatches[g], linewidth=2.0, alpha=0.75)
ax[0].set_ylabel("Mismatches (#)")

ax[1].plot(x_ax, PROBS, color="red", linewidth=3.0)
ax[1].plot(x_ax, PROBS_UW, color="green", linewidth=3.0)
ax[1].set_ylabel("Log(P)")
ax[1].set_title("Path Likelihood by Iteration")

ax[2].plot(x_ax, WEIGHTS)
ax[2].set_ylabel("Observations (#)")
ax[2].set_title("Observations Removed by Iteration")

# Add lines
for k, v in best_by_id.items():
    ax[0].axvline(v["p_"], color='k', linewidth=2.0, alpha=v["id"]/100.0)
    ax[1].axvline(v["p_"], color='k', linewidth=2.0, alpha=v["id"]/100.0)
    ax[2].axvline(v["p_"], color='k', linewidth=2.0, alpha=v["id"]/100.0)
    #ax[2].axvline(r, color='k', linestyle='--')
plt.show()

"""
#TODO None for ARGS
if HITS and REFS:
    con_mat, snp_mat, seen_mat, master_mat, full_confusion = gretel.confusion_matrix(PATHS, VCF_h, HITS, REFS, REF_NAMES, VCF_h["N"], ARGS.master)

    if not ARGS.quiet:
        print "#\tname\tsites\trate\tbestit0\trefd\tlogl\tmap"
    RECOVERIES = []
    RECOV_PCTS = []

    for i in range(0, len(con_mat)):
        recovered = 0.0
        at = None
        for j in range(0, len(con_mat[i])):
            if con_mat[i][j] > recovered:
                recovered = con_mat[i][j]
                at = j
        if at != None:
            RECOVERIES.append((at, recovered))
            if not ARGS.quiet:
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
        RECOV_PCTS.append(recovered)
    if ARGS.quiet:
        print "%.2f %.2f %.2f" % (min(RECOV_PCTS), max(RECOV_PCTS), np.mean(RECOV_PCTS))


        # Add likelihood
        ax[1].plot(x_ax, PATH_PROBS, color="red", linewidth=3.0)
        ax[1].plot(x_ax, PATH_PROBS_UW, color="green", linewidth=3.0)
        ax[1].set_ylabel("Log(P)")
        ax[1].set_title("Path Likelihood by Iteration")

        #ax[2].plot(x_ax, PATH_PROBS_RATIO, color="blue", linewidth=3.0)
        ax[1].set_xlabel("Iteration (#)")

        # Add recoveries
        for r in RECOVERIES:
            ax[0].axvline(r[0], color='k', linewidth=2.0, alpha=r[1]/100.0)
            ax[1].axvline(r[0], color='k', linewidth=2.0, alpha=r[1]/100.0)
            #ax[2].axvline(r, color='k', linestyle='--')
        plt.show()

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

