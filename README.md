# q2-confluence: a [rachis](https://rachis.org) plugin for use with QIIME 2

> ⚠️ This plugin is in early testing (alpha release) phase. Interfaces will almost certainly change. Be skeptical of results - testing is still in progress.

`q2-confluence` (registered name `confluence`; `qiime confluence --help`) integrates 16S rRNA amplicon and metagenome taxonomic profiles in the feature space of one GTDB release's genomes and compares the studies on the samples both profiled.
The metagenome profile comes from the companion plugin [mosh-sylph](https://github.com/gregcaporaso/mosh-sylph).

The documentation is the tutorial, [thermo-hec-tutorial](https://gregcaporaso.github.io/thermo-hec-tutorial/) ([source](https://github.com/gregcaporaso/thermo-hec-tutorial)), which also says how to create an environment with both plugins; `qiime confluence <action> --help` is the action reference.

Development: create the tutorial's environment, then `make dev`, `make test` and `make lint` from this directory.

The code, tests and usage examples were written by Claude and reviewed function by function by @gregcaporaso; the commit history says what was requested and what the assistant implemented.
Questions and bug reports go to the [QIIME 2 Forum](https://forum.qiime2.org).
