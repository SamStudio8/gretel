Protocol
========

**Gretel** provides a command line tool for the recovery of haplotypes.
We recommend the following protocol.

Read Alignment
--------------

**Gretel** requires your reads to be aligned to a common reference. This is to
ensure that reads share a co-ordinate system, on which we can call for variants
and recover haplotypes. The reference itself is of little consequence, though
dropped reads will lead to evidence to be unavailable to Gretel.

Construction of a *de novo* consensus assembly for a metagenome is left as an exercise
for the reader. Align the reads to your assembly (`bowtie2`, `minimap2` etc.).
Sort and index the alignment BAM.

Variant Calling
---------------

**Gretel** is robust to sequencing error and misalignment noise, thus the
calling of variants need not be carefully conducted. Typically we have used `samtools`,
but for our own Gretel pipeline, we have aggressively called all heterogenous sites
in an alignment as a SNP using the `snpper` tool in our `gretel-test repository
<https://github.com/SamStudio8/gretel-test>`_.

For somewhat questionable reasoning, we currently require a compressed and indexed VCF: ::

    bgzip <my.vcf>
    tabix <my.vcf.gz>

Invocation of Gretel
--------------------
As described in the README, Gretel is invoked as follows: ::

    gretel <my.sort.bam> <my.vcf.gz> <contig> [-s 1startpos] [-e 1endpos] [--master master.fa] [-o output_dir]

You must provide your sorted BAM, compressed VCF, and the name of the contig on which
to recover haplotypes. Use `-s` and `-e` to specify the positions on the aligned reads between which
to recover haplotypes from your metagenome.

By default, Gretel will output a FASTA containing the recovered SNPs, in order, for each haplotype.
Providing an optional "master" FASTA sequence will permit Gretel to "fill in" the non-SNP positions
(*i.e.* the positions between `-s` and `-e` that do not appear in the VCF) with the nucleotide from
the pseudo-reference.

Gretel Outputs
--------------

out.fasta
~~~~~~~~~
A **FASTA** containing each of the recovered sequences, in the order they were found.
Each sequence is named `<iteration>__-<log10 likelihood>`. Sequences are not wrapped.

gretel.crumbs
~~~~~~~~~~~~~

Additionally, Gretel outputs a whimsically named *crumbs* file, containing some potentially
interesting metadata, as well as a record of each recovered haplotype.
The first row is a comment containing the following (in order):

* The number of SNPs across the region of interest
* The number of 'crumbs': paired observations added to the Hansel matrix
* The number of 'slices': reads with at least one observation added to the Hansel matrix
* The chosen value of `L` for the `L`'th order Markov chain

The rest of the file contains tab-delimited metadata for each recovered haplotype:

* The iteration number, starting from 0
* The number of times this haplotype was returned
* The *weighted* likelihood of the haplotype, given the Hansel matrix at the time the haplotype was recovered (comma-sep for each time the haplotype was returned)
* The *unweighted* likelihood of the haplotype, given the Hansel matrix at the time the reads were parsed (comma-sep for each time the haplotype was returned)
* The haplotype magnitude: total number of observations removed from the Hansel matrix by the reweighting mechanism

In practice, we rank with the **weighted** likelihoods to discern the haplotypes most likely to exist in the metagenome.
One may attempt to use the *unweighted* likelihoods as a means to compare the abundance, or read support, **between the returned haplotypes** (*i.e.* not necessarily the metagenome as a whole).
