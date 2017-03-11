import sys
from math import log,log10,exp
import random

import vcf
import pysam
import numpy as np

from hansel import Hansel
import util

#TODO Should the denom of the conditional use the unique variants at i-l or i?
#TODO Util to parse known input and return SNP seq

def reweight_hansel_from_path(hansel, path, ratio):
    #TODO Reweight only the L pairs?
    """
    Given a completed path, reweight the applicable pairwise observations in the Hansel structure.

    Parameters
    ----------
    hansel : :py:class:`hansel.hansel.Hansel`
        The Hansel structure currently being explored by Gretel.

    path : list{str}
        The ordered sequence of selected variants.

    ratio : float
        The proportion of evidence to remove from each paired observation that
        was considered to recover the provided path.

        It is recommended this be the smallest marginal distribution observed across selected variants.

        *i.e.* For each selected variant in the path, note the value of the
        marginal distribution for the probability of observing that particular
        variant at that genomic position. Parameterise the minimum value of
        those marginals.

    Returns
    -------
    Spent Observations : float
        The sum of removed observations from the Hansel structure.
    """

    size = 0
    for i in range(0, len(path)):
        for j in range(0, i+1+1):
            # Reduce read supports
            if i >= len(path)-1:
                size += hansel.reweight_observation(path[i], path[j], i, i+1, ratio)
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
                size += hansel.reweight_observation(path[t_i], path[t_j], t_i, t_j, ratio)
    sys.stderr.write("[RWGT] Ratio %.3f, Removed %.1f\n" % (ratio, size))
    return size


## INPUT OUTPUT ###############################################################
def process_vcf(vcf_path, contig_name, start_pos, end_pos):
    """
    Parse a VCF to extract the genomic positions of called variants.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.

    contig_name : str
        Name of the target contig on which variants were called.

    start_pos : int
        The 1-indexed genomic position from which to begin considering variants.

    end_pos : int
        The 1-indexed genomic position at which to stop considering variants.

    Returns
    -------
    Gretel Metastructure : dict
        A collection of structures used for the execution of Gretel.
        The currently used keys are:
            N : int
                The number of observed SNPs
            snp_fwd : dict{int, int}
                A reverse lookup from the n'th variant, to its genomic position on the contig
            snp_rev : dict{int, int}
                A forward lookup to translate the n'th genomic position to its i'th SNP rank
            region : list{int}
                A masked representation of the target contig, positive values are variant positions
    """

    # Open the VCF
    vcf_records = vcf.Reader(open(vcf_path))
    n_snps = 0
    snp_reverse = {}
    snp_forward = {}
    region = np.zeros(end_pos + 1, dtype=int)
    i = 0
    for record in vcf_records.fetch(contig_name, 0, end_pos):
        if record.POS < start_pos:
            continue
        if record.POS > end_pos:
            continue

        n_snps += 1
        region[record.POS] = 1
        snp_reverse[i] = record.POS
        snp_forward[record.POS] = i
        i += 1

    return {
        "N": n_snps,
        "snp_fwd": snp_forward,
        "snp_rev": snp_reverse,
        "region": region,
    }

def process_bam(vcf_handler, bam_path, contig_name, start_pos, end_pos, L, use_end_sentinels, n_threads):
    """
    Initialise a Hansel structure and load variants from a BAM.

    Parameters
    ----------
    vcf_handler : dict{str, any}
        Variant metadata, as provided by :py:func:`gretel.gretel.process_vcf`.

    bam_path : str
        Path to the alignment BAM.

    contig_name : str
        The name of the contig for which to recover haplotypes.

    start_pos : int
        The 1-indexed genomic position from which to begin considering variants.

    end_pos : int
        The 1-indexed genomic position at which to stop considering variants.

    L : int
        The Gretel `L-parameter`, controlling the number of positions back
        from the head of the current path (including the head) to consider
        when calculating conditional probabilities.

    use_end_sentinels : boolean
        Whether or not to append an additional pairwise observation between
        the final variant on a read towards a sentinel.

    n_threads : int
        Number of threads to spawn for reading the BAM

    Returns
    -------
    Gretel Metastructure : dict
        A collection of structures used for the execution of Gretel.
        The currently used keys are:
            read_support : :py:class:`hansel.hansel.Hansel`
                The Hansel structure.
            read_support_o : :py:class:`hansel.hansel.Hansel`
                A copy of the Hansel structure stored with the intention of not reweighting its observations.
            meta : dict{str, any}
                A dictionary of metadata returned from the BAM parsing, such as
                a list of the number of variants that each read spans.
    """

    #NOTE(samstudio8)
    # Could we optimise for lower triangle by collapsing one of the dimensions
    # such that Z[m][n][i][j] == Z[m][n][i + ((j-1)*(j))/2]


    meta = util.load_from_bam(bam_path, contig_name, start_pos, end_pos, vcf_handler, use_end_sentinels, n_threads)
    hansel = Hansel(meta["hansel"], ['A', 'C', 'G', 'T', 'N', "-", "_"], ['N', "_"], L=L)

    if hansel.L == 0:
        hansel.L = meta["L"]
        sys.stderr.write("[NOTE] Setting Gretel.L to %d\n" % hansel.L)

    return {
        "read_support": hansel,
        "read_support_o": hansel.copy(),
        "meta": meta,
    }

## PATH GENERATION ############################################################

def append_path(path, next_m, next_v):
    """
    Append a selected variant to a given path.
    .. deprecated:: 1.0
        This method is somewhat of a stub.
        It is likely to be deprecated at no notice in future.

    Parameters
    ----------
    path : list{str}
        The current sequence of variants representing a path (haplotype) in progress.

    next_m : str
        The symbol to append to the path.

    next_v : float
        The marginal probability of `next_m` at the current position.

    Raises
    ------
    Exception
        Raised if `next_m` is None.

    """
    #TODO(samstudio8) This is somewhat of a pointless stub, now.
    #TODO(samstudio8) Probably a bit gross as it has side effects on path...
    #TODO(samstudio8) Could probably raise an Exception for any next_m not in hansel.symbols?
    if next_m is not None:
        path.append(next_m)
    else:
        raise Exception("Cowardly refusing to append None as a nucleotide. Cheerio.")

def generate_path(n_snps, hansel, original_hansel):
    """
    Explore and generate the most likely path (haplotype) through the observed Hansel structure.

    Parameters
    ----------
    n_snps : int
        The number of variants.

    hansel : :py:class:`hansel.hansel.Hansel`
        The Hansel structure currently being explored by Gretel.

    original_hansel : :py:class:`hansel.hansel.Hansel`
        A copy of the Hansel structure created by Gretel, before any reweighting.

    Returns
    -------
    Path : list{str} or None
        The sequence of variants that represent the completed path (or haplotype), or None
        if one could not be successfully constructed.

    Path Probabilities : dict{str, float}
        The `unweighted` (orignal Hansel) and `weighted` (current Hansel) joint
        probabilities of the variants in the returned path occurring together
        in the given order.

    Minimum Marginal : float
        The smallest marginal distribution observed across selected variants.
    """

    # Cross the metahaplome in a greedily, naive fashion to establish a base path
    # This seeds the rest of the path generation (we might want to just select
    #   a random path here in future)

    running_prob = 0.0
    running_prob_uw = 0.0
    current_path = ['_'] # start with the dummy
    marginals = []

    # Find path
    sys.stderr.write("*** ESTABLISH ***\n")
    for snp in range(1, n_snps+1):
        #sys.stderr.write("\t*** ***\n")
        #sys.stderr.write("\t[SNP_] SNP %d\n" % snp)

        # Get marginal and calculate branch probabilities for each available
        # mallele, given the current path seen so far
        # Select the next branch and append it to the path
        curr_branches = hansel.get_edge_weights_at(snp, current_path)
        #sys.stderr.write("\t[TREE] %s\n" % curr_branches)
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

        selected_edge_weight = hansel.get_marginal_of_at(next_m, snp)
        marginals.append(selected_edge_weight) #NOTE This isn't a log, as it is used as a ratio later

        running_prob += log10(selected_edge_weight)
        running_prob_uw += log10(original_hansel.get_marginal_of_at(next_m, snp))
        append_path(current_path, next_m, next_v)

    return current_path, {"unweighted": running_prob_uw, "weighted": running_prob}, min(marginals)

