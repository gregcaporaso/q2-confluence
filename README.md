# q2-confluence: a [rachis](https://rachis.org) plugin for use with QIIME 2

⚠️ This plugin is in early testing (alpha release) phase. Interfaces will almost certainly change. Be skeptical of results - testing is still in progress.

`q2-confluence` integrates 16S rRNA amplicon and metagenome taxonomic profiles of the same samples in the feature space of one GTDB release's genomes.
It is designed to be a convenience layer over stock QIIME 2 [Actions](https://news.rachis.org/en/latest/glossary#term-action).
The metagenome profiles can come from the companion plugin [mosh-sylph](https://github.com/gregcaporaso/mosh-sylph).

See the [thermo-hec-tutorial](https://gregcaporaso.github.io/thermo-hec-tutorial/) ([source](https://github.com/gregcaporaso/thermo-hec-tutorial)) for documentation, including installation and usage instructions.
Run `qiime confluence --help` after installation to learn what actions are available.

Questions and bug reports should be directed to the [QIIME 2 Forum](https://forum.qiime2.org).

## AI disclosure

The code, tests and usage examples were written by Claude and reviewed function by function by @gregcaporaso; the commit history says what was requested and what the assistant implemented.
As of this writing (28 October 2026) all code has been reviewed, while usage example and unit test review is still pending.
Functions that are unreviewed are tagged in the [source code](https://github.com/gregcaporaso/q2-confluence) with the comment line:

```
# TODO: this is AI-generated code; human review pending
```

[This search](https://github.com/search?q=repo%3Agregcaporaso%2Fq2-confluence+%23+TODO%3A+this+is+AI-generated+code%3B+human+review+pending&type=code) will direct you to all of those lines.
