"""Given a BAM, a contig, and an ending genomic position, aggressively call for
variants and generate a placeholder VCF.

Thanks to @linsalrob for the initial argparsification of gretel-snpper.
"""
import sys

import numpy as np
import pysam
import argparse

def main():
    parser = argparse.ArgumentParser('Aggressively call for variants and generate a VCF', epilog='NOTE: Coordinates are 1-based as they are for samtools')
    parser.add_argument('--bam', help='bam of reads aligned to (psuedo)-reference', required=True)
    parser.add_argument('--contig', help='name of contig to generate a VCF for', required=True)
    parser.add_argument('-s', help='start (default = 1)', type=int, default=1)
    parser.add_argument('-e', help='end (default = length of the reference)', type=int)
    parser.add_argument('--depth', help='number of reads that must feature a read to call that base as a possible variant (default = 0)', type=int, default=0)
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)

    if not args.e:
        args.e = bam.lengths[bam.references.index(args.contig)]

    # convert 1-indexed numbers to 0 indexed numbers for pysam
    args.s = args.s - 1

    counts = np.array(bam.count_coverage(contig=args.contig, start=args.s, stop=args.e, quality_threshold=0, read_callback='nofilter'))

    COUNT_SENSITIVITY = args.depth

    vcf_h = [
        "##fileformat=VCFv4.2",
    ]
    vcf = []

    sites = (counts > COUNT_SENSITIVITY).sum(axis=0)
    for i, s in enumerate(sites):
        if s > 1:
            vcf.append([
                args.contig,
                i+1+args.s,
                '.',
                'A',
                'C,T,G',
                0,
                '.',
                "INFO"
            ])


    for r in vcf_h:
        print(r)
    for r in vcf:
        print("\t".join([str(s) for s in r]))
