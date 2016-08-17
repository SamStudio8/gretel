gretel
======
A (work-in-progress) algorithm for recovering haplotypes from metagenomes.
Sister to `Hansel<https://github.com/SamStudio8/hansel>`_.

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

Citation
--------

Please cite us so we can continue to make useful software! ::

    Nicholls, S.M., Aubrey, W., de Grave, K., Schietgat, L., Creevey, C.J. and Clare, A., 2016. Advances in the recovery of haplotypes from the metagenome. bioRxiv, p.067215.

::

    @article {Nicholls067215,
        author = {Nicholls, Samuel M and Aubrey, Wayne and de Grave, Kurt and Schietgat, Leander and Creevey, Chris J and Clare, Amanda},
        title = {Advances in the recovery of haplotypes from the metagenome},
        year = {2016},
        doi = {10.1101/067215},
        publisher = {Cold Spring Harbor Labs Journals},
        URL = {http://biorxiv.org/content/early/2016/08/02/067215},
        eprint = {http://biorxiv.org/content/early/2016/08/02/067215.full.pdf},
        journal = {bioRxiv}
    }

License
-------
Hansel and Gretel are distributed under the MIT license, see LICENSE.
