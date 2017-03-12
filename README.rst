Gretel
======
An algorithm for recovering haplotypes from metagenomes.
Sister to `Hansel
<https://github.com/SamStudio8/hansel>`_.

What is it?
-----------

**Gretel** is a Python package providing a command line tool for the recovery of haplotypes
from metagenomic data sets. **Gretel** parses an alignment of reads into a **Hansel** matrix
and uses the evidence of SNP pairs observed to appear on the same reads to probabilistically
reconstruct the most likely haplotypes.

**Gretel** uses an L'th order Markov chain model to reconstruct likely sequences
of variants that constitute haplotypes in the real metagenome.
Our approach involves graph-like traversal of the data within the **Hansel** matrix.
Edges are probabilitically weighted based on the evidence on the reads, as well as
the haplotype as it has been reconstructed so far.

What can I use it for?
----------------------

**Gretel** is designed to recover haplotypes from your data set, without the need for
setting (or optimisation) of any parameters.
**Gretel** does not require a priori knowledge of your input data (such as its contents, or
the true number of haplotypes) and makes no assumptions
regarding the distributions of alleles at variant sites and uses the available evidence
from the aligned reads without altering or discarding the observed varations.

Why should I use it?
--------------------

**Gretel** is the first tool capable of recovering haplotypes from metagenomes.
Whilst tools exist for analogous haplotyping problems, such as the assembly of
viral quasispecies, typically these tools rely on overlap approaches that create
too many unranked haplotypes. **Gretel** is capable of ranking the haplotypes it
outputs by their likelihood.

**Gretel** requires no parameters and our approach is robust to sequencing error
and misalignment noise.

Requirements
------------
::

    $ pip install numpy hanselx pysam PyVCF

Install
-------
::

    $ pip install gretel

Usage
-----
You will require a sorted BAM containing your reads, aligned to some pseudo-reference.
You can use any sequence as your reference, such as a consensus assembly of the
metagenomic reads, or a known strain reference (such as HIV-1).
You must bgzip and tabix your VCF.

::

    $ gretel <bam> <vcf.gz> <contig> -s <1-start> -e <1-end> --master <master.fa> -o <outdir>


Citation
--------
Paper pending...

License
-------
Hansel and Gretel are distributed under the MIT license, see LICENSE.
