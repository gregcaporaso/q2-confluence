(tutorial)=
# Integrating 16S and shotgun metagenome profiles: the tutorial

<!-- Draft by Claude, pending @gregcaporaso's rewrite (no prose of theirs yet).
Every command is a usage example rendered for each interface; the data is the
thermo-hec tutorial subset (ten compost samples with both 16S V4 amplicon and
shotgun metagenome data, GTDB r232). How the data and the reference artifacts
were prepared is on the [data preparation](data-preparation.md) page.
Drive IDs are those of the 2026-10-16 upload (tutorial-data/README.md in the parent project). -->

This tutorial integrates a 16S rRNA amplicon profile and a shotgun metagenome profile of the same samples on one GTDB release, with the `confluence` plugin and stock QIIME 2 actions.
It covers what you run on your own data: classifying and filtering the amplicon features, profiling the metagenome reads, integrating both at one taxonomic level, and comparing the two profile types on the samples they share.
[How the tutorial data and the reference artifacts were prepared](data-preparation.md) is a separate page.

<!-- @gregcaporaso, reminder for the prose: say that the abundance-weighted results on this 10% subset equal the full-depth run (bacon job 308 versus jobs 300 to 307: weighted UniFrac pseudo-F 2.61 against 2.62, Bray-Curtis 3.21 against 3.20) while the presence-based metrics and the feature counts move (unweighted UniFrac 4.62 against 3.83, Jaccard 4.93 against 4.58, 2,020 against 2,385 features) because sylph detects fewer genomes at this depth. -->

<!-- @gregcaporaso, reminder for the prose: explain why a reader's numbers differ slightly from the ones shown. The PERMANOVA p-values are permutation estimates and the random forest draws its folds and trees afresh each run, so p-values wobble in the third decimal and the accuracy by a sample or two (0.90 and 0.95 on two runs of this page on the same inputs); the pseudo-F values and the feature summary are deterministic. No random seed is set on purpose: copied commands would reuse it. -->

## Working directory

Create a directory for the tutorial and change into it; every file below is written there.

```shell
mkdir confluence-tutorial
cd confluence-tutorial
```

## Sample metadata

The sample metadata of the ten tutorial samples: bin location, sample type, composting time point, core temperature in °C, and the sample UUID.

:::{describe-usage}
:scope: confluence-tutorial

sample_metadata = use.init_metadata_from_url(
    'sample-metadata',
    'https://drive.usercontent.google.com/download?export=download&confirm=t&id=1atS6Xt5npYZv0YBLrWiXVGMb4YO0u0RT')
:::

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(input=sample_metadata),
    use.UsageOutputNames(visualization='sample-metadata'))
:::

The integration actions take the study metadata as an `ImmutableMetadata` artifact, so the same table is also provided imported.

:::{describe-usage}
study_metadata = use.init_artifact_from_url(
    'study-metadata',
    'https://drive.usercontent.google.com/download?export=download&confirm=t&id=1PJQOCPLgOJ_Jeb-_ZmEN1enFzqLhbH5S')
:::

## Reference data

Both profile types are placed on GTDB r232.
The GTDB bacterial tree collapsed to one tip per genus was prepared once for the release ([how](data-preparation.md#reference)); the amplicon classifiers come from their publishers, and the sylph database from the sylph authors.

:::{describe-usage}
genus_tree = use.init_artifact_from_url(
    'gtdb-r232-genus-tree',
    'https://drive.usercontent.google.com/download?export=download&confirm=t&id=1sgwhDnhSXIbliFb-Fdi3o1u2NExqbovD')

gtdb_classifier = use.init_artifact_from_url(
    'gtdb-r232-classifier',
    'https://zenodo.org/records/21619532/files/gtdb-r232.0-2026.7.0-classifier.qza?download=1')

silva_classifier = use.init_artifact_from_url(
    'silva-144-v4-classifier',
    'https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/V4-515f-806r/uniform/SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza')
:::

The sylph database of the GTDB r232 species representatives is about 20 GB; `fetch-sylph-database` downloads it with the matching sylph-tax taxonomy.

:::{describe-usage}
sylph_database, = use.action(
    use.UsageAction(plugin_id='confluence', action_id='fetch_sylph_database'),
    use.UsageInputs(release='232'),
    use.UsageOutputNames(database='sylph-gtdb-r232'))
:::

:::{tip} An artifact cache for the large artifacts
:class: dropdown
Every action that reads a `.qza` file unzips it first, and the sylph database is 20 GB.
An [artifact cache](https://use.qiime2.org/en/stable/tutorials/use-the-artifact-cache/) keeps artifacts unzipped under a key, referenced as `cache:key` in place of a file path, so the database is written once and never zipped or unzipped again.
On the command line, the alternative to the command above is:

```shell
qiime tools cache-create --cache cache

qiime confluence fetch-sylph-database \
  --p-release 232 \
  --use-cache cache \
  --o-database cache:sylph-gtdb-r232
```

Then pass `--i-database cache:sylph-gtdb-r232 --use-cache cache` to `profile-reads` and `query-reads` below.
The tutorial's own data stays in `.qza` files, with its provenance.
:::

## Amplicon data

The amplicon side starts from the denoised ASV table and sequences of the ten samples.

:::{describe-usage}
asv_table = use.init_artifact_from_url(
    'asv-table',
    'https://drive.usercontent.google.com/download?export=download&confirm=t&id=1-HH4892h693WAQ37YvrBjV_3P7xA3UgF')

asv_sequences = use.init_artifact_from_url(
    'asv-sequences',
    'https://drive.usercontent.google.com/download?export=download&confirm=t&id=1LzOioUEnzP4nbOsV9v75VGaqfgseagbi')
:::

The GTDB classifier gives every ASV a GTDB taxonomy; the SILVA classifier is used only to find organelle sequences (mitochondria, chloroplasts), which GTDB cannot name.

:::{describe-usage}
asv_taxonomy_gtdb, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_sklearn'),
    use.UsageInputs(classifier=gtdb_classifier, reads=asv_sequences, n_jobs=8),
    use.UsageOutputNames(classification='asv-taxonomy-gtdb'))

asv_taxonomy_silva, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_sklearn'),
    use.UsageInputs(classifier=silva_classifier, reads=asv_sequences, n_jobs=8),
    use.UsageOutputNames(classification='asv-taxonomy-silva'))
:::

Organelles are removed with the SILVA taxonomy, then the table is restricted to Bacteria with the GTDB taxonomy, because the reference tree is the bacterial GTDB tree and `integrate` refuses features that are not tips of it.

:::{describe-usage}
asv_table_noorganelles, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(table=asv_table, taxonomy=asv_taxonomy_silva,
                    exclude='f__Mitochondria,o__Chloroplast,d__Eukaryota',
                    mode='contains'),
    use.UsageOutputNames(filtered_table='asv-table-noorganelles'))

asv_table_bacteria, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(table=asv_table_noorganelles, taxonomy=asv_taxonomy_gtdb,
                    include='d__Bacteria', mode='contains'),
    use.UsageOutputNames(filtered_table='asv-table-bacteria'))
:::

## Metagenome data

The metagenome side starts from paired-end reads that are already trimmed, free of control samples and human reads, and subsampled to 10% ([how](data-preparation.md#metagenome)).

:::{describe-usage}
reads = use.init_artifact_from_url(
    'wgs-subset-reads',
    'https://drive.usercontent.google.com/download?export=download&confirm=t&id=1fkzsYoWAIeGOJDQrdqGcf5k5rNfoFR7E')
:::

sylph sketches the reads once per sample; the sketches serve every profile or query against any database.

:::{describe-usage}
sketches, = use.action(
    use.UsageAction(plugin_id='confluence', action_id='sketch_reads'),
    use.UsageInputs(reads=reads, threads=16),
    use.UsageOutputNames(sketches='sketches'))
:::

:::{tip} Sketches in the cache
:class: dropdown
The sketches are a few GB; the cache form of the step above is:

```shell
qiime confluence sketch-reads \
  --i-reads wgs-subset-reads.qza \
  --p-threads 16 \
  --use-cache cache \
  --o-sketches cache:sketches
```
:::

`profile-reads` runs sylph profile and sylph-tax taxprof: per sample, the share of each detected genome corrected for genome size (the taxonomic abundance, the measure to integrate), the share of reads (the sequence abundance), the genome-keyed taxonomy, and per-sample statistics (genomes and genera detected, median ANI, unknown share).

:::{describe-usage}
sylph_taxonomic_abundance, sylph_sequence_abundance, sylph_taxonomy, sylph_stats = use.action(
    use.UsageAction(plugin_id='confluence', action_id='profile_reads'),
    use.UsageInputs(sketches=sketches, database=sylph_database, threads=16),
    use.UsageOutputNames(taxonomic_abundance='sylph-taxonomic-abundance',
                         sequence_abundance='sylph-sequence-abundance',
                         taxonomy='sylph-taxonomy', stats='sylph-stats'))
:::

:::{tip} Profiling from the cache
:class: dropdown
With the database and the sketches in the cache:

```shell
qiime confluence profile-reads \
  --i-sketches cache:sketches \
  --i-database cache:sylph-gtdb-r232 \
  --p-threads 16 \
  --use-cache cache \
  --o-taxonomic-abundance sylph-taxonomic-abundance.qza \
  --o-sequence-abundance sylph-sequence-abundance.qza \
  --o-taxonomy sylph-taxonomy.qza \
  --o-stats sylph-stats.qza
```
:::

The per-sample statistics are metadata and tabulate like any other.

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(input=use.view_as_metadata('sylph-stats-md', sylph_stats)),
    use.UsageOutputNames(visualization='sylph-stats'))
:::

As on the amplicon side, the profile is restricted to Bacteria.

:::{describe-usage}
sylph_bacteria, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(table=sylph_taxonomic_abundance, taxonomy=sylph_taxonomy,
                    include='d__Bacteria', mode='contains'),
    use.UsageOutputNames(filtered_table='sylph-bacteria'))
:::

## Integration

The two profiles enter as collections keyed by a label that names the profile type in every output.
`summarize-inputs` shows which sample IDs the profiles and the metadata share before anything is merged.

:::{describe-usage}
tables = use.construct_artifact_collection(
    'tables', {'amplicon': asv_table_bacteria, 'metagenome': sylph_bacteria})
taxonomies = use.construct_artifact_collection(
    'taxonomies', {'amplicon': asv_taxonomy_gtdb, 'metagenome': sylph_taxonomy})
metadata = use.construct_artifact_collection(
    'metadata', {'study': study_metadata})

use.action(
    use.UsageAction(plugin_id='confluence', action_id='summarize_inputs'),
    use.UsageInputs(tables=tables, metadata=metadata),
    use.UsageOutputNames(visualization='inputs-summary'))
:::

`integrate` collapses each profile to the chosen level (6 is genus with the seven GTDB ranks; 5 is family), drops the features without a name at that level, renormalizes each profile to proportions, and merges them.
Merged samples are named `<sample-id>-<key>`; the integrated metadata records `original-sample-id`, `profile-type` and `unassigned-share` beside the study metadata.

:::{describe-usage}
relative_table, taxonomy, integrated_metadata = use.action(
    use.UsageAction(plugin_id='confluence', action_id='integrate'),
    use.UsageInputs(tables=tables, taxonomies=taxonomies, tree=genus_tree,
                    metadata=metadata, level=6),
    use.UsageOutputNames(relative_table='genus-relative-table',
                         taxonomy='genus-taxonomy',
                         integrated_metadata='genus-metadata'))
:::

## Paired-sample report

On the samples both profile types cover, `paired-report` computes Bray-Curtis, Jaccard, unweighted and weighted UniFrac on the proportions as they are (no rarefaction), with PERMANOVA by profile type and a PCoA scatter plot for each; a random forest classifying the profile type with nested cross-validation; and a per-feature summary, tabulated over every feature and also returned as metadata.

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='confluence', action_id='paired_report'),
    use.UsageInputs(table=relative_table, tree=genus_tree,
                    metadata=integrated_metadata,
                    permutations=999, cv=10, n_estimators=500, n_jobs=8),
    use.UsageOutputNames(report='genus-paired-report',
                         distance_matrices='genus-distance-matrices',
                         feature_importance='genus-feature-importance',
                         feature_summary='genus-feature-summary'))
:::

For the family-level view, repeat the integration and the report with `level=5` and the family tree, prepared like the genus tree: `https://drive.usercontent.google.com/download?export=download&confirm=t&id=1MM2gsJFqSF-Ec7nHcqvKVGReliEjDM_R`.
