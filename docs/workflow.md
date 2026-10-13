(workflow)=
# The integration workflow, command by command

<!-- Draft by Claude, pending @gregcaporaso's rewrite (2026-10-13: no prose of
theirs yet). This page is the command sequence of the thermo-hec baseline (ten
compost samples with both 16S V4 amplicon and shotgun metagenome data, GTDB
r232), with one-line captions. Data download links are placeholders until the
tutorial subset is published. -->

Stock QIIME 2 and MOSHPIT actions prepare each side; `confluence` actions start at the reference tree.
Commands assume a working directory with `16S/`, `wgs/` and `reference/` subdirectories and the `q2-confluence` environment active.

## Reference artifacts (once per GTDB release)

The bacterial tree and the genome-keyed taxonomy of a release, collapsed to the level of the comparison.

```shell
qiime tools import \
  --type 'Phylogeny[Rooted]' \
  --input-path bac120_r232.tree \
  --output-path reference/gtdb-r232-bac120-tree.qza

gunzip -k bac120_taxonomy_r232.tsv.gz
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path bac120_taxonomy_r232.tsv \
  --output-path reference/gtdb-r232-bac120-taxonomy.qza

qiime confluence gtdb-tree-collapse \
  --i-tree reference/gtdb-r232-bac120-tree.qza \
  --i-taxonomy reference/gtdb-r232-bac120-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-tree reference/gtdb-r232-genus-tree.qza
```

Level 6 is genus with the seven GTDB ranks; 5 is family and 7 species.
The same collapsed tree serves every table collapsed to that level.

## Amplicon side (stock actions)

Starting from the denoised ASV table and sequences.
The GTDB classifier is the rachis release for r232; the SILVA classifier is the 2026.7 V4 uniform classifier, used only to find organelles, which GTDB cannot name.

```shell
qiime feature-classifier classify-sklearn \
  --i-classifier reference/gtdb-r232.0-2026.7.0-classifier.qza \
  --i-reads 16S/asv-sequences.qza \
  --p-n-jobs 8 \
  --o-classification 16S/asv-taxonomy-gtdb.qza

qiime feature-classifier classify-sklearn \
  --i-classifier reference/SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza \
  --i-reads 16S/asv-sequences.qza \
  --p-n-jobs 8 \
  --o-classification 16S/asv-taxonomy-silva.qza

qiime taxa filter-table \
  --i-table 16S/asv-table.qza \
  --i-taxonomy 16S/asv-taxonomy-silva.qza \
  --p-exclude f__Mitochondria,o__Chloroplast,d__Eukaryota \
  --p-mode contains \
  --o-filtered-table 16S/asv-table-noorganelles.qza
```

On the thermo-hec data the organelle filter removed 329 of 35,731 ASVs and 0.5% of the reads.

## Metagenome side

Import, trim, drop the controls, remove human reads: stock actions.

```shell
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path wgs-manifest.tsv \
  --output-path wgs/demux.qza

mosh fastp process-seqs \
  --i-sequences wgs/demux.qza \
  --p-length-required 90 --p-cut-mean-quality 30 --p-cut-tail \
  --o-processed-sequences wgs/reads.qza \
  --o-reports wgs/fastp-reports.qza

mosh demux filter-samples \
  --i-demux wgs/reads.qza \
  --m-metadata-file controls.tsv \
  --p-exclude-ids \
  --o-filtered-demux wgs/reads-nocontrols.qza

mosh quality-control filter-reads-human-pangenome \
  --i-reads wgs/reads-nocontrols.qza \
  --p-threads 16 \
  --o-filtered-reads wgs/reads-nohost.qza \
  --o-reference-index reference/human-pangenome-bowtie2.qza
```

The index build downloads GRCh38 and the HPRC pangenome and takes hours the first time; pass `--i-index reference/human-pangenome-bowtie2.qza` on later runs.

Then the sylph profile, two `confluence` actions.

```shell
qiime confluence fetch-sylph-database \
  --p-release 232 \
  --o-database reference/sylph-gtdb-r232.qza

qiime confluence profile-reads \
  --i-reads wgs/reads-nohost.qza \
  --i-database reference/sylph-gtdb-r232.qza \
  --p-threads 16 \
  --o-taxonomic-abundance wgs/sylph-taxonomic-abundance.qza \
  --o-sequence-abundance wgs/sylph-sequence-abundance.qza \
  --o-taxonomy wgs/sylph-taxonomy.qza \
  --o-stats wgs/sylph-stats.qza
```

The taxonomic abundance (genome-size corrected) is the measure to integrate; the sequence abundance (read share) is kept alongside.
`sylph-stats` holds, per sample, the genomes and genera detected, the median ANI and the unknown share.

## Integration

Collections are passed as directories whose file names are the keys.

```shell
mkdir -p tables taxonomies metadata
cp 16S/asv-table-noorganelles.qza tables/amplicon.qza
cp wgs/sylph-taxonomic-abundance.qza tables/metagenome.qza
cp 16S/asv-taxonomy-gtdb.qza     taxonomies/amplicon.qza
cp wgs/sylph-taxonomy.qza        taxonomies/metagenome.qza

qiime tools import \
  --type ImmutableMetadata \
  --input-path sample-metadata.tsv \
  --output-path metadata/study.qza

qiime confluence summarize-inputs \
  --i-tables tables/ \
  --i-metadata metadata/ \
  --o-visualization inputs-summary.qzv

qiime confluence integrate \
  --i-tables tables/ \
  --i-taxonomies taxonomies/ \
  --i-tree reference/gtdb-r232-genus-tree.qza \
  --i-metadata metadata/ \
  --p-level 6 \
  --o-relative-table genus/relative-table.qza \
  --o-taxonomy genus/taxonomy.qza \
  --o-integrated-metadata genus/metadata.qza
```

`summarize-inputs` shows which sample IDs the sources share; `integrate` pairs samples by identical IDs and names the merged samples `<id>-amplicon` and `<id>-metagenome`, with `original-sample-id`, `source` and `unassigned-share` in the integrated metadata.
Features without a name at the level are dropped from every source, and each source is renormalized to proportions, so counts (the amplicon table) and proportions (the sylph table) integrate alike.

## Paired-sample report

```shell
qiime confluence paired-report \
  --i-table genus/relative-table.qza \
  --i-tree reference/gtdb-r232-genus-tree.qza \
  --i-metadata genus/metadata.qza \
  --p-permutations 999 \
  --p-cv 10 \
  --p-n-estimators 500 \
  --p-n-jobs 8 \
  --o-report genus/paired-report.qzv \
  --o-distance-matrices genus/distance-matrices/ \
  --o-feature-importance genus/feature-importance.qza
```

Bray-Curtis, Jaccard, unweighted and weighted UniFrac are computed on the proportions as they are, without rarefaction.
The report holds, per metric, PERMANOVA by source and a PCoA scatter plot colored by source, then the random forest's confusion matrix and the per-feature summary.
Repeat the collapse, integration and report with `--p-level 5` for the family-level view.
