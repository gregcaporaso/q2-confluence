# q2-confluence

A [QIIME 2](https://qiime2.org) plugin by Greg Caporaso that integrates 16S rRNA amplicon and shotgun metagenome taxonomic profiles on one GTDB release.
The registered plugin name is `confluence` (`qiime confluence --help`); the Python package is `q2_confluence`.

Status: in development, not yet released; the repository is private for now.
Actions marked deprecated at registration are provisional and name the plugin they are expected to move to.

## Actions

The workflow, in order (`qiime confluence <action> --help` for the details):

| Action | Role | Status |
|---|---|---|
| `fetch-sylph-database` | Pre-built sylph database of the GTDB species representatives (226 or 232) with its sylph-tax taxonomy | provisional, intended for MOSHPIT |
| `profile-reads` | sylph sketch, profile and taxprof on paired-end reads: genome-keyed taxonomic and sequence abundance tables, taxonomy, per-sample stats | provisional, intended for MOSHPIT |
| `gtdb-tree-collapse` | A reference tree collapsed to one tip per taxon at a level, named like `taxa collapse` IDs; once per release and level | provisional, intended for q2-phylogeny |
| `summarize-inputs` | Which sample IDs the tables and metadata of the studies share | stays |
| `merge-metadata` | Union or intersection of several studies' metadata | stays |
| `relative-frequency` | Proportions from counts or from a relative table whose features were filtered; `integrate` applies it to every source | provisional, intended for q2-feature-table |
| `integrate` | Every source (counts or proportions) collapsed to the level, unnamed features dropped, renormalized, merged with `original-sample-id` and `source` metadata columns; sample IDs become `<id>-<key>` | stays |
| `paired-report` | On the samples every source profiled: Bray-Curtis, Jaccard and UniFrac on the proportions (no rarefaction), PERMANOVA by source and a PCoA scatter plot per metric, random forest with nested cross-validation, per-feature summary, one report | stays |
| `summarize-paired-features` | The per-feature summary on its own | stays |

The amplicon side (denoising, GTDB classification, the SILVA-based organelle filter) and the metagenome read preparation (trimming, host read removal) are stock QIIME 2 and MOSHPIT actions and are described in the documentation rather than wrapped here.
Every action carries a usage example that the test suite executes.

## Development installation

Create the environment (QIIME 2 2026.7 `qiime2` distribution, sylph, and the documentation toolchain), then install the plugin into it:

```shell
conda env create -n q2-confluence-dev --file environment-files/q2-confluence-qiime2-dev.yml
conda activate q2-confluence-dev
make dev
qiime dev refresh-cache
qiime confluence --help
```

On Apple Silicon, prefix the `conda env create` with `CONDA_SUBDIR=osx-64` and run `conda config --env --set subdir osx-64` after activating.

Run the tests and the linters with:

```shell
make test
make lint
```

## Documentation

The documentation is a MyST book built with q2doc, in the layout of the [gut-to-soil tutorial](https://github.com/qiime2/gut-to-soil-tutorial), and is published on ReadTheDocs from `.readthedocs.yaml`.
Build it locally, in the environment above with the plugin installed, with:

```shell
q2doc refresh-cache
make html
make serve
```

## About

`q2-confluence` was [created from the QIIME 2 plugin template](https://develop.qiime2.org/en/stable/plugins/tutorials/create-from-template.html) on 2026-10-09.
It is a community plugin: questions and bug reports go to the [QIIME 2 Forum](https://forum.qiime2.org), where the developers will see them; see the [QIIME 2 Library](https://library.qiime2.org) for the community-plugin support policy.
