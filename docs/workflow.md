(workflow)=
# The integration workflow, command by command

<!-- Draft by Claude, pending @gregcaporaso's rewrite (no prose of theirs yet).
This page is the command sequence of the thermo-hec baseline (ten compost
samples with both 16S V4 amplicon and shotgun metagenome data, GTDB r232),
with one-line captions. -->

Stock QIIME 2 and MOSHPIT actions prepare each side; `confluence` actions start at the reference tree.
Commands assume the `q2-confluence` environment is active.

## Working directory

Everything the tutorial downloads or produces, the artifact cache included, lives under one directory.

```shell
mkdir confluence-tutorial
cd confluence-tutorial
mkdir 16S wgs reference
```

## Data

The tutorial subset: the ten compost samples with both assays.
The shotgun reads are already trimmed, control-free and human-read-free, and subsampled to 10% (2.4 to 4.8 million pairs per sample), a depth at which the sylph profile is indistinguishable from full depth by weighted and unweighted UniFrac.
The 16S artifacts are the unfiltered DADA2 outputs of the whole study (199 samples); the organelle filter below is part of the workflow.

```shell
D='https://drive.usercontent.google.com/download?export=download&confirm=t&id='
wget -O wgs/wgs-subset-reads.qza "${D}1fkzsYoWAIeGOJDQrdqGcf5k5rNfoFR7E"
wget -O 16S/asv-table.qza "${D}1-mcm1_ngVVuysqgJyURKaOTEI0OV1vee"
wget -O 16S/asv-sequences.qza "${D}1uOcdxD5WhiH83qVNnvq55kkwMCKHd58B"
wget -O sample-metadata.tsv "${D}1BgqmKoAfQpUmlKA9Mc86FK8uCnB7fX-l"
md5sum wgs/wgs-subset-reads.qza 16S/asv-table.qza 16S/asv-sequences.qza sample-metadata.tsv
```

| File | Size | md5 | Artifact UUID (`qiime tools peek`) |
|---|---|---|---|
| `wgs/wgs-subset-reads.qza` | 4.8 GB | `884d47f6948d010356cde49594238b3f` | `09901ee0-9bca-4825-b94a-20c953d97cef` |
| `16S/asv-table.qza` | 2.7 MB | `48d9f13d671e30714b3aad3f8665ec98` | `c24a9747-46fe-456e-bb74-6bc5799a7710` |
| `16S/asv-sequences.qza` | 2.7 MB | `cdd603a2e8dae3d9f2fb643a81ae50b2` | `b74df7ea-a4a8-4aa4-95c6-869f955d4444` |
| `sample-metadata.tsv` | 37 KB | `f0523491b674cb19383ebbc3f7780778` | |

The reads artifact carries the provenance of the import, trimming, control removal, host-read removal and subsampling described under "Metagenome side".

Two of the artifacts in this workflow are large (the sylph database is about 20 GB, the read sketches a few GB), and every action that reads a `.qza` unzips it first.
An [artifact cache](https://use.qiime2.org/en/stable/tutorials/use-the-artifact-cache/) keeps artifacts unzipped under a key, referenced as `cache:key` in place of a file path, so the large ones are written once and never unzipped again.
The commands below keep the database and the sketches in a cache called `cache` inside the tutorial directory and pass `--use-cache cache` to the actions that read them; the downloaded artifacts stay as `.qza` files, with their provenance.
For your own data, `qiime tools cache-import` imports a large input (a read manifest, a database directory) straight into a cache key, so it is never zipped at all.

```shell
qiime tools cache-create --cache cache
```

The reference files come from their sources: the GTDB r232 bacterial tree and genome taxonomy, the GTDB r232 naive Bayes classifier from the rachis Zenodo record, and the SILVA 144 V4 uniform classifier from SILVA.

```shell
wget -P reference https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_r232.tree
wget -P reference https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv.gz
wget -O reference/gtdb-r232.0-2026.7.0-classifier.qza 'https://zenodo.org/records/21619532/files/gtdb-r232.0-2026.7.0-classifier.qza?download=1'
wget -P reference https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/V4-515f-806r/uniform/SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza
```

## Reference artifacts (once per GTDB release)

The bacterial tree and the genome-keyed taxonomy of a release, collapsed to the level of the comparison.

```shell
qiime tools import \
  --type 'Phylogeny[Rooted]' \
  --input-path reference/bac120_r232.tree \
  --output-path reference/gtdb-r232-bac120-tree.qza

gunzip -k reference/bac120_taxonomy_r232.tsv.gz
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path reference/bac120_taxonomy_r232.tsv \
  --output-path reference/gtdb-r232-bac120-taxonomy.qza

qiime confluence gtdb-tree-collapse \
  --i-tree reference/gtdb-r232-bac120-tree.qza \
  --i-taxonomy reference/gtdb-r232-bac120-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-tree reference/gtdb-r232-genus-tree.qza
```

Level 6 is genus with the seven GTDB ranks; 5 is family and 7 species.
The same collapsed tree serves every table collapsed to that level.
On r232 the collapse takes about ten seconds per level: 6,048 families, 34,834 genera, 189,801 species.

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

How the subset was prepared, with stock actions, from the study's raw reads: import, trim, drop the controls, remove human reads.
These steps are recorded in the provenance of `wgs/wgs-subset-reads.qza`; a reader working on the subset starts at the sylph profile below.

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

The index build downloads GRCh38 and the HPRC pangenome and takes hours the first time (on a 16-core server: 4 h to build at 50 GB peak, then 5 h 42 min to filter 393 million pairs at 58 GB peak, removing 3.8% of them); pass `--i-index reference/human-pangenome-bowtie2.qza` on later runs.
The subset was then drawn with `qiime demux subsample-paired --p-fraction 0.1`.

Then the sylph profile, three `confluence` actions, on the subset (minutes; the full data took 3 min at 4.7 GB peak).

```shell
qiime confluence fetch-sylph-database \
  --p-release 232 \
  --use-cache cache \
  --o-database cache:sylph-gtdb-r232

qiime confluence sketch-reads \
  --i-reads wgs/wgs-subset-reads.qza \
  --p-threads 16 \
  --use-cache cache \
  --o-sketches cache:sketches

qiime confluence profile-reads \
  --i-sketches cache:sketches \
  --i-database cache:sylph-gtdb-r232 \
  --p-threads 16 \
  --use-cache cache \
  --o-taxonomic-abundance wgs/sylph-taxonomic-abundance.qza \
  --o-sequence-abundance wgs/sylph-sequence-abundance.qza \
  --o-taxonomy wgs/sylph-taxonomy.qza \
  --o-stats wgs/sylph-stats.qza
```

The sketches are made once per sample and serve every later profile or query against any database.

The taxonomic abundance (genome-size corrected) is the measure to integrate; the sequence abundance (read share) is kept alongside.
`sylph-stats` holds, per sample, the genomes and genera detected, the median ANI and the unknown share.

## Integration

The reference tree is the bacterial GTDB tree, so each source is first restricted to Bacteria; `integrate` refuses any feature that is not a tip of the tree, and both the amplicon and the sylph profile hold archaeal genera.
Collections are passed as directories, one artifact per key.

```shell
mkdir -p tables taxonomies metadata
qiime taxa filter-table \
  --i-table 16S/asv-table-noorganelles.qza \
  --i-taxonomy 16S/asv-taxonomy-gtdb.qza \
  --p-include d__Bacteria \
  --p-mode contains \
  --o-filtered-table tables/amplicon.qza
qiime taxa filter-table \
  --i-table wgs/sylph-taxonomic-abundance.qza \
  --i-taxonomy wgs/sylph-taxonomy.qza \
  --p-include d__Bacteria \
  --p-mode contains \
  --o-filtered-table tables/metagenome.qza
cp 16S/asv-taxonomy-gtdb.qza taxonomies/amplicon.qza
cp wgs/sylph-taxonomy.qza    taxonomies/metagenome.qza

qiime tools import \
  --type ImmutableMetadata \
  --input-path sample-metadata.tsv \
  --output-path metadata/study.qza

printf 'amplicon\nmetagenome\n' | tee tables/.order > taxonomies/.order
printf 'study\n' > metadata/.order
```

The `.order` file in each directory names the keys; without one, QIIME 2 keys the collection by file name, and the sources would be called `amplicon.qza` and `metagenome.qza` in every output.

```shell

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
