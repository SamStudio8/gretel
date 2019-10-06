import sys
from math import log,log10,exp
import random

import numpy as np

from hansel import Hansel
from . import util

#TODO Should the denom of the conditional use the unique variants at i-l or i?
#TODO Util to parse known input and return SNP seq

def reweight_hansel_from_path(hansel, path, ratio):
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

    """
    # Old re-implementation sans flip
    for i in range(0, len(path)-1):
        for j in range(0, i+1+1):
            # Reduce read supports
            if i == j:
                continue
            size += hansel.reweight_observation(path[i], path[j], i, j, ratio)
    return size

    # Reduce adjacent evidence pairs
    for i in range(len(path)-1):
        size += hansel.reweight_observation(path[i], path[i+1], i, i+1, ratio)

    # Reduce other evidence pairs
    for j in range(1, len(path)):
        for i in range(0, j-1):
            size += hansel.reweight_observation(path[i], path[j], i, j, ratio)

    # Reduce other non-evidence pairs
    # I have no idea why this works so well, so we'll need to have a think about it
    # before we put it in Gretel proper...
    #for j in range(1, len(path)):
    #    for i in range(0, j-1):
    #        size += hansel.reweight_observation(path[j], path[i], j, i, ratio)
    #        pass

    # Reweight the rest of the matrix because we can at least explain that
    hansel.reweight_matrix( ratio / (hansel.L/10) )

    sys.stderr.write("[RWGT] Ratio %.3f, Removed %.1f\n" % (ratio, size))
    return size
    """

    # Let's keep the RW system as-is for now...
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

## PATH GENERATION ############################################################

def generate_path(n_snps, hansel, original_hansel, debug_hpos=None):
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
        The `hp_original` (orignal Hansel) and `hp_current` (current Hansel) joint
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
    current_path = [ hansel.symbols_d['_'] ] # start with the dummy
    marginals = []

    # Find path
    sys.stderr.write("[NOTE] *Establishing next path\n")
    for snp in range(1, n_snps+1):
        #sys.stderr.write("\t*** ***\n")
        #sys.stderr.write("\t[SNP_] SNP %d\n" % snp)

        dh_flag = False
        if debug_hpos:
            if snp in debug_hpos:
                dh_flag = True

        # Get marginal and calculate branch probabilities for each available
        # mallele, given the current path seen so far
        # Select the next branch and append it to the path
        curr_branches = hansel.get_edge_weights_at(snp, current_path, debug=dh_flag)
        #sys.stderr.write("\t[TREE] %s\n" % curr_branches)
        # Return the symbol and probability of the next base to add to the
        # current path based on the best marginal
        next_v = 0.0
        next_m = None

        if debug_hpos:
            if snp in debug_hpos:
                print(curr_branches)

        for symbol in curr_branches:
            if str(symbol) == "total":
                continue
            if next_m is None:
                next_v = curr_branches[symbol]
                next_m = symbol
            elif curr_branches[symbol] > next_v:
                next_v = curr_branches[symbol]
                next_m = symbol

        if next_m == None:
            sys.stderr.write('''[NOTE] Unable to select next branch from SNP %d to %d
       By design, Gretel will attempt to recover haplotypes until a hole in the graph has been found.
       Recovery will intentionally terminate now.\n''' % (snp-1, snp))
            return None, None, None

        selected_edge_weight = hansel.get_marginal_of_at(next_m, snp)
        marginals.append(selected_edge_weight) #NOTE This isn't a log, as it is used as a ratio later

        running_prob += log10(selected_edge_weight)
        running_prob_uw += log10(original_hansel.get_marginal_of_at(next_m, snp))
        current_path.append(next_m)

    return current_path, {"hp_original": running_prob_uw, "hp_current": running_prob}, min(marginals)

