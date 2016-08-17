gretel
======
A (work-in-progress) algorithm for recovering potential haplotypes from metagenomes

Housekeeping
------------

Dependencies
~~~~~~~~~~~~
::

    $ pip install numpy hanselx pysam PyVCF matplotlib

Install
~~~~~~~
::

    $ pip install gretel

virtualenv
~~~~~~~~~~

`matplotlib` misbehaves somewhat when you attempt to install it in a fashion that
isolates it from the system packages, to get around this one can cheat and include
the system packages in the `virtualenv`.:: 

    $ virtualenv gretel --system-site-packages


Usage
-----
::

    $ gretel <bam> <vcf.gz> <contig> -e <1-end> --master <master.fa>


To plot some graphs given a known FASTA of haplotypes and the Gretel crumbs: ::

    $ makeblastdb -in out.fasta -dbtype nucl -out out.blast
    $ blastn -query known_haplotypes.fa -db out.blast -outfmt 6 | sort -k2 -n > out.hit
    $ gretel-crumbs <out.hit> <gretel.crumbs>

License
-------
Hansel and Gretel are distributed under the MIT license, see LICENSE.
