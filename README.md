# q2-confluence

A [QIIME 2](https://qiime2.org) plugin by Greg Caporaso that integrates 16S rRNA amplicon and shotgun metagenome taxonomic profiles on one GTDB release.
The registered plugin name is `confluence` (`qiime confluence --help`); the Python package is `q2_confluence`.

Status: in development, not yet released; the repository is private for now.
Actions marked deprecated at registration are provisional and name the plugin they are expected to move to.

## Actions

The workflow, in order (`qiime confluence <action> --help` for the details):

| Action | Role | Status |
|---|---|---|
| `gtdb-tree-collapse` | A reference tree collapsed to one tip per taxon at a level, named like `taxa collapse` IDs, plus a log (`ImmutableMetadata` keyed by tip: original node name, genome count, genome IDs); once per release and level, for stock diversity actions on tables collapsed with `taxa collapse` (the workflow itself integrates at the genome level and does not collapse) | provisional, intended for q2-phylogeny |
| `summarize-inputs` | Which sample IDs the tables and metadata of the studies share | stays |
| `relative-frequency` | Proportions from counts or from a relative table whose features were filtered; `integrate` applies it to every source | provisional, intended for q2-feature-table |
| `integrate` | Every source (counts or proportions, keyed by GTDB genome ID) filtered to the tips of the GTDB tree (`phylogeny filter-table`, counts in the metadata and a warning past a tenth lost), renormalized and merged in the feature space of the reference genomes, with `original-sample-id`, `study-id`, `features-before`, `features-after` and `features-kept-share` metadata columns and the union of the studies' metadata (`metadata merge`); sample IDs become `<id>-<key>`; the taxonomy is GTDB's own file, used directly downstream | stays |
| `paired-report` | On the samples every source profiled: Bray-Curtis, Jaccard and UniFrac on the proportions (no rarefaction), PERMANOVA by study-id and a PCoA scatter plot per metric, random forest with nested cross-validation, per-feature summary (every feature, searchable, with the taxonomy beside each ID when given), one report plus the summary as metadata | stays |
| `summarize-paired-features` | The per-feature summary (importance, rank, per-study mean and prevalence, log2 ratio, direction) as `ImmutableMetadata`; `paired-report` tabulates it | stays |

The amplicon side (denoising, the SILVA-based organelle filter, closed-reference clustering of the ASVs against the GTDB SSU reference, which keys them by genome) and the metagenome read preparation (trimming, host read removal) are stock QIIME 2 and MOSHPIT actions, and the metagenome profile comes from the `sylph` plugin ([mosh-sylph](https://github.com/gregcaporaso/mosh-sylph), thin wrappers around sylph; its `profile-to-tables` keys the tables by GTDB genome ID). The tutorial uses both plugins in one environment.
Every action carries a usage example that the test suite executes.

## Formats

| Format | Role | Status |
|---|---|---|
| `GTDBNewickFormat` | A GTDB reference tree as released (`bac120_r232.tree`, `ar53_r232.tree`): drops the root label UniFrac cannot parse and quotes the tips so every reader keeps their underscores; `qiime tools import --type 'Phylogeny[Rooted]' --input-format GTDBNewickFormat` takes the file unchanged, and the transformer to the stock Newick view drops the quoted root label without a branch length (`'100:d__Bacteria'`) that UniFrac's parser (unifrac-binaries 1.5, run by `diversity beta-phylogenetic`) aborts on; nothing else in the tree changes | provisional, intended for q2-types or q2-phylogeny |

Known limitations of the genome feature space, documented in the tutorial: only about half of GTDB r232's bacterial species representatives have a 16S in the SpeciesReps SSU reference, so an ASV whose genome has none can only land on a relative's genome at a lower identity threshold; and `vsearch cluster-features-closed-reference` keeps vsearch's first hit at or above the threshold (q2-vsearch does not expose `maxaccepts`), so among genomes with identical V4 regions the assignment is arbitrary.

## Development installation

The plugin targets the unreleased 2026.10 epoch, so the base environment is the `qiime2` distribution's development environment, created as in the [QIIME 2 development environment guide](https://develop.qiime2.org/en/stable/plugins/how-to-guides/set-up-development-environment.html) from the distribution's `dev/latest/passed` environment file for your platform; the plugin's extras (sylph, sylph-tax, the documentation toolchain, the linters) are layered on with `conda env update`.

Linux (x86_64):

```shell
conda env create \
  -n q2-confluence-dev-$(date "+%Y-%m-%d") \
  --file https://raw.githubusercontent.com/qiime2/distributions/dev/latest/passed/rachis-qiime2-linux-64-conda.yml
```

macOS on Apple Silicon (the distribution has no arm64 build; the osx-64 build runs under Rosetta):

```shell
CONDA_SUBDIR=osx-64 conda env create \
  -n q2-confluence-dev-$(date "+%Y-%m-%d") \
  --file https://raw.githubusercontent.com/qiime2/distributions/dev/latest/passed/rachis-qiime2-osx-64-conda.yml

conda activate q2-confluence-dev-$(date "+%Y-%m-%d")

conda config \
  --env \
  --set subdir osx-64
```

macOS on Intel uses the same osx-64 file without `CONDA_SUBDIR`.
Then, in the activated environment:

```shell
conda env update \
  -n q2-confluence-dev-$(date "+%Y-%m-%d") \
  --file environment-files/q2-confluence-dev-extras.yml

make dev
qiime dev refresh-cache
qiime confluence --help
```

The dev environment files are rebuilt as the distribution's packages change; recreate the environment every few weeks, as the guide suggests.

Run the tests and the linters with:

```shell
make test
make lint
```

## Documentation

The tutorial that uses this plugin, with the data-preparation record, is its own repository, [thermo-hec-tutorial](https://github.com/gregcaporaso/thermo-hec-tutorial).
This repository carries no documentation site; the action reference is `qiime confluence <action> --help` and the usage example each action registers.

## About

`q2-confluence` was [created from the QIIME 2 plugin template](https://develop.qiime2.org/en/stable/plugins/tutorials/create-from-template.html) on 2026-10-09.
It is a community plugin: questions and bug reports go to the [QIIME 2 Forum](https://forum.qiime2.org), where the developers will see them; see the [QIIME 2 Library](https://library.qiime2.org) for the community-plugin support policy.
