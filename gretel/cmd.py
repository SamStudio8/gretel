import argparse
import numpy as np
import sys
import os

from . import gretel
from . import util

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

    parser.add_argument("--master", default=None, help="Master sequence (will be used to fill in homogeneous gaps in haplotypes, otherwise --gapchar)") #TODO Use something other than N? Should probably be a valid IUPAC
    parser.add_argument("--gapchar", default="N", help="Character to fill homogeneous gaps in haplotypes if no --master [default N]")
    parser.add_argument("--delchar", default="", help="Character to output in haplotype for deletion (eg. -) [default is blank]")

    parser.add_argument("--quiet", default=False, action='store_true', help="Don't output anything other than a single summary line.")
    #parser.add_argument("--sentinels", default=False, action='store_true', help="Add additional sentinels for read ends [default:False][EXPERIMENTAL]")
    parser.add_argument("-o", "--out", default=".", help="Output directory [default .]")
    parser.add_argument("-@", "--threads", type=int, default=1, help="Number of BAM iterators [default 1]")

    parser.add_argument("--debugreads", type=str, default="", help="A newline delimited list of read names to output debug data when parsing the BAM")
    parser.add_argument("--debugpos", type=str, default="", help="A newline delimited list of 1-indexed genomic positions to output debug data when parsing the BAM")
    parser.add_argument("--debughpos", type=str, default=",", help="A comma delimited list of 1-indexed SNP positions to output debug data when predicting haplotypes")

    parser.add_argument("--dumpmatrix", type=str, default=None, help="Location to dump the Hansel matrix to disk")
    parser.add_argument("--dumpsnps", type=str, default=None, help="Location to dump the SNP positions to disk")

    ARGS = parser.parse_args()

    debug_hpos = []
    if ARGS.debughpos:
        for x in ARGS.debughpos.split(","):
            try:
                debug_hpos.append( int(x) )
            except:
                pass

    if ARGS.end == -1:
        ARGS.end = util.get_ref_len_from_bam(ARGS.bam, ARGS.contig)
        sys.stderr.write("[NOTE] Setting end_pos to %d" % ARGS.end)

    debug_reads = set([])
    if ARGS.debugreads:
        debug_fofn = open(ARGS.debugreads)
        for line in debug_fofn:
            debug_reads.add(line.strip())

    debug_pos = set([])
    if ARGS.debugpos:
        debug_fofn = open(ARGS.debugpos)
        for line in debug_fofn:
            debug_pos.add(int(line.strip()))

    VCF_h = gretel.process_vcf(ARGS.vcf, ARGS.contig, ARGS.start, ARGS.end)
    if ARGS.dumpsnps:
        snp_fh = open(ARGS.dumpsnps, 'w')
        for k in sorted(VCF_h["snp_fwd"].keys()):
            snp_fh.write("%d\t%d\t%d\n" % (VCF_h["snp_fwd"][k]+1, k, k-ARGS.start+1))
        snp_fh.close()

    # Could we optimise for lower triangle by collapsing one of the dimensions
    # such that Z[m][n][i][j] == Z[m][n][i + ((j-1)*(j))/2]
    hansel = util.load_from_bam(ARGS.bam, ARGS.contig, ARGS.start, ARGS.end, VCF_h, n_threads=ARGS.threads, debug_reads=debug_reads, debug_pos=debug_pos)
    original_hansel = hansel.copy()

    if ARGS.dumpmatrix:
        hansel.save_hansel_dump(ARGS.dumpmatrix)

    # Check if there is a gap in the matrix
    # TODO(samstudio8) Ideally we would do this IN the bam worker threads such
    #                  that a thread could raise an error and halt the whole
    #                  process if we already know we'll yield a matrix we cannot
    #                  actually work with, but this will do for now...
    for i in range(0, VCF_h["N"]+1):
        marginal = hansel.get_counts_at(i)

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

    PATHS = {}

    # Spew out exciting information about the SNPs
    if not ARGS.quiet:
        print ("i\tpos\tgap\tA\tC\tG\tT\tN\t-\t_\ttot")
        last_rev = 0
        for i in range(0, VCF_h["N"]+1):
            marginal = hansel.get_counts_at(i)
            marginal = {str(x): marginal[x] for x in marginal}
            snp_rev = 0
            if i > 0:
                snp_rev = VCF_h["snp_rev"][i-1]
            print ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (
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
            ))
            last_rev = snp_rev


    # Make some genes
    SPINS = ARGS.paths
    ongoing_mag = 0
    for i in range(0, SPINS):
        init_path, init_prob, init_min = gretel.generate_path(VCF_h["N"], hansel, original_hansel, debug_hpos=debug_hpos)
        if init_path == None:
            break
        current_path = init_path

        MIN_REMOVE = 0.01 # 1%
        if init_min < MIN_REMOVE:
            sys.stderr.write("[RWGT] Ratio %.10f too small, adjusting to %.3f\n" % (init_min, MIN_REMOVE))
            init_min = MIN_REMOVE
        rw_magnitude = gretel.reweight_hansel_from_path(hansel, init_path, init_min)

        #TODO Horribly inefficient.
        current_path_str = "".join([str(x) for x in current_path])
        if current_path_str not in PATHS:
            PATHS[current_path_str] = {
                "hp_current": [],
                "hp_original": [],
                "i": [],
                "i_0": i,
                "n": 0,
                "magnitude": 0,
                "hansel_path": current_path,
            }
        PATHS[current_path_str]["n"] += 1
        PATHS[current_path_str]["i"].append(i)
        PATHS[current_path_str]["magnitude"] += rw_magnitude
        PATHS[current_path_str]["hp_current"].append(init_prob["hp_current"])
        PATHS[current_path_str]["hp_original"].append(init_prob["hp_original"])

    # Make some pretty pictures
    dirn = ARGS.out + "/"
    fasta_out_fh = open(dirn+"out.fasta", "w")
    hfasta_out_fh = open(dirn+"snp.fasta", "w")

    if ARGS.master:
        master_fa = util.load_fasta(ARGS.master)
        master_seq = master_fa.fetch(master_fa.references[0])
    else:
        master_seq = [' '] * ARGS.end

    for p in sorted(PATHS, key=lambda x: PATHS[x]["i_0"]):
        p = PATHS[p]
        path = p["hansel_path"]
        i = p["i_0"]

        seq = list(master_seq[:])
        for j, mallele in enumerate(path[1:]):
            snp_pos_on_master = VCF_h["snp_rev"][j]
            try:
                if mallele == hansel.symbols_d["-"]:
                    # It's a deletion, don't print a SNP
                    seq[snp_pos_on_master-1] = ARGS.delchar
                else:
                    seq[snp_pos_on_master-1] = mallele
            except IndexError:
                print (path, len(seq), snp_pos_on_master-1)
                sys.exit(1)

        # Coerce HanselSymbols to str
        to_write = "".join(str(x) for x in seq[ARGS.start-1 : ARGS.end])
        if not ARGS.master:
            to_write = to_write.replace(' ', ARGS.gapchar)

        fasta_out_fh.write(">%d__%.2f\n" % (i, p["hp_current"][0])) #TODO hp_current or hp_original?
        fasta_out_fh.write("%s\n" % to_write)

        hfasta_out_fh.write(">%d__%.2f\n" % (i, p["hp_current"][0])) #TODO hp_current or hp_original?
        hfasta_out_fh.write("%s\n" % "".join([str(x) for x in path[1:]]))
    fasta_out_fh.close()
    hfasta_out_fh.close()

    #TODO datetime, n_obs, n_slices, avg_obs_len, L, n_paths, n_avg_loglik
    crumb_file = open(dirn+"gretel.crumbs", "w")
    crumb_file.write("# %d\t%d\t%d\t%.2f\n" % (
        VCF_h["N"],
        hansel.n_crumbs,
        hansel.n_slices,
        hansel.L,
    ))

    for p in sorted(PATHS, key=lambda x: PATHS[x]["hp_current"][0], reverse=True):
        p = PATHS[p]
        crumb_file.write("%d\t%d\t%s\t%s\t%.2f\n" % (
                p["i_0"],
                p["n"],
                ",".join(["%.2f" % x for x in p["hp_current"]]),
                ",".join(["%.2f" % x for x in p["hp_original"]]),
                p["magnitude"],
        ))
