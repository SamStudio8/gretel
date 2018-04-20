History
=======

0.0.90
------
Resolves a bug whereby SNPs are incorrectly parsed from the BAM if either:
* its quality score is below 13
* the read is overlapped by its primary mate
Well covered data sets need not be overly affected by the additional noise that
may have been introduced, but the problem is more noticeable with low coverage
and you may wish to reapply Gretel to affected data. Sorry.


0.0.81
------
* Add warning and advice when an entry in Hansel is missing evidence.
* Make the 'Unable to select' warning sound much less bad because it is normal.

0.0.8
-----
* Docs
* Deprecate `gretel-crumbs` command

0.0.7
-----
* Further improvements to parallel read processing
* Add `-` symbol to enable support for deletions

0.0.6b
------
* Fix setting of `L` parameter

0.0.6
-----
* MULTIPROCESSING
* Re-write read handling, again

0.0.5
-----
* `-s` and `-e` introduced to allow specification of positions between which
  to recover haplotypes
* Attempt some basic indel handling
* Fix a bug where the master sequence was altered by the output of each
  reported haplotype

0.0.4
-----
* Add experimental `--sentinels` option
* Improve docs

0.0.3
-----
* Hansel is now seperate from Gretel
* [Hansel] `get_marginal_at` is `now get_counts_at`
* [Hansel] `selext_next_edge_at` deprecated
* Gene recovery and likelihood plots are now on seperate panels
* Re-write methods to add observations to matrix to be less awful to read
* Drop `--hit` and `--gene` options to verification
* Replace verification script to `gretel-crumbs` command

0.0.2
-----
* Improve documentation.
* Provide `util` subpackage for filling `Hansel` structure with BAM observations.
* Explicitly provide possible symbols to `Hansel`.
* Improve plotting
* Remove `process_hits` and `process_refs` as these are no longer needed.
* Rename `establish_path` to `generate_path`
* Rename `add_ignore_support3` to `reweight_hansel_from_graph` so we have some sort of indication of what it does.
* Altered Sphinx configuation.

0.0.1
-----
* Import repository from `claw`.
