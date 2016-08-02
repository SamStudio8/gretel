gretel
======
A (work-in-progress) algorithm for recovering potential haplotypes from metagenomes

Housekeeping
------------

virtualenv
~~~~~~~~~~

`matplotlib` misbehaves somewhat when you attempt to install it in a fashion that
isolates it from the system packages, to get around this one can cheat and include
the system packages in the `virtualenv`.:: 

    $ virtualenv gretel --system-site-packages

Dependencies
~~~~~~~~~~~~
::

    $ pip install numpy hanselx pysam PyVCF matplotlib

Install
~~~~~~~
::

    $ pip install gretel

Usage
-----
::

    $ gretel <bam> <vcf.gz> <contig> -e <1-end> --master <master.fa>
