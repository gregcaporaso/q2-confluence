(data-preparation)=
# How the tutorial data and the reference artifacts were prepared

<!-- Draft by Claude, pending @gregcaporaso's rewrite. These steps are the
record of what produced the hosted artifacts; a reader of the tutorial does
not run them. The subsetting of the amplicon data to the ten tutorial samples
is not shown: it is in the artifacts' provenance. -->

None of this is part of the tutorial: it records how the artifacts the [tutorial](tutorial.md) downloads were produced, so that the same can be done for another GTDB release or another data set.
Every step is also in the provenance of the artifacts themselves.

(reference)=
## Reference

The GTDB r232 bacterial tree and genome taxonomy, from GTDB, imported and collapsed to one tip per genus (level 6 of the seven GTDB ranks) and per family (level 5) with `gtdb-tree-collapse`.
The collapsed trees are the hosted reference artifacts; the tutorial downloads the genus tree.
On r232 each collapse takes about ten seconds: 34,834 genera, 6,048 families.

```shell
wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_r232.tree
wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv.gz
gunzip -k bac120_taxonomy_r232.tsv.gz

qiime tools import \
  --type 'Phylogeny[Rooted]' \
  --input-path bac120_r232.tree \
  --output-path gtdb-r232-bac120-tree.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path bac120_taxonomy_r232.tsv \
  --output-path gtdb-r232-bac120-taxonomy.qza

qiime confluence gtdb-tree-collapse \
  --i-tree gtdb-r232-bac120-tree.qza \
  --i-taxonomy gtdb-r232-bac120-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-tree gtdb-r232-genus-tree.qza

qiime confluence gtdb-tree-collapse \
  --i-tree gtdb-r232-bac120-tree.qza \
  --i-taxonomy gtdb-r232-bac120-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-tree gtdb-r232-family-tree.qza
```

The amplicon classifiers and the sylph database are not re-hosted: the tutorial fetches the GTDB r232 naive Bayes classifier from the rachis Zenodo record, the SILVA 144 V4 uniform classifier from SILVA, and the sylph database and taxonomy from the sylph authors with the `sylph` plugin's `fetch-database` and `tax-download`.

(amplicon)=
## Amplicon data

The 16S rRNA V4 amplicons (515F/806R) of the study were sequenced on an Illumina run of 199 samples and denoised with DADA2 in the QIIME 2 amplicon distribution, giving the ASV table and the ASV sequences.
The tutorial artifacts are those outputs restricted to the ten samples with a metagenome profile (`feature-table filter-samples` and `filter-seqs`; the step is in their provenance).
Classification, the organelle filter and the Bacteria filter are part of the tutorial, because they are what a reader runs on their own data.

(metagenome)=
## Metagenome data

The shotgun reads of the same ten samples were imported from a fastq manifest, trimmed with fastp, cleared of the run's control samples, cleared of human reads with the MOSHPIT human pangenome filter, and subsampled to 10%.
The subset holds 2.4 to 4.8 million read pairs per sample (median 4.1 million); at that depth the sylph profile is indistinguishable from the full depth by weighted and unweighted UniFrac.
The hosted reads artifact carries the provenance of every step.

```shell
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path wgs-manifest.tsv \
  --output-path demux.qza

mosh fastp process-seqs \
  --i-sequences demux.qza \
  --p-length-required 90 --p-cut-mean-quality 30 --p-cut-tail \
  --o-processed-sequences reads.qza \
  --o-reports fastp-reports.qza

mosh demux filter-samples \
  --i-demux reads.qza \
  --m-metadata-file controls.tsv \
  --p-exclude-ids \
  --o-filtered-demux reads-nocontrols.qza

mosh quality-control filter-reads-human-pangenome \
  --i-reads reads-nocontrols.qza \
  --p-threads 16 \
  --o-filtered-reads reads-nohost.qza \
  --o-reference-index human-pangenome-bowtie2.qza

qiime demux subsample-paired \
  --i-sequences reads-nohost.qza \
  --p-fraction 0.1 \
  --o-subsampled-sequences wgs-subset-reads.qza
```

The human pangenome index downloads GRCh38 and the HPRC pangenome and takes hours the first time (on a 16-core server: 4 h to build at 50 GB peak, then 5 h 42 min to filter 393 million pairs at 58 GB peak, removing 3.8% of them); pass `--i-index human-pangenome-bowtie2.qza` on later runs.

(metadata)=
## Sample metadata

The study's sample sheet was reduced to the ten tutorial samples and to five columns, renamed in kebab-case: `bin-location`, `sample-type`, `composting-time-point`, `core-temp-c` (the compost core temperature in °C; two samples, the fresh excrement and the soil, have none) and `sample-uuid`.
The same table is hosted as a TSV and imported as an `ImmutableMetadata` artifact for the integration actions.
