# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Bool, Choices, Citations, Collection, Int, Metadata,
                           Plugin, Range, Str, Visualization)
from q2_types.feature_data import FeatureData, Importance, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.metadata import ImmutableMetadata
from q2_types.per_sample_sequences import PairedEndSequencesWithQuality
from q2_sample_classifier import ClassifierPredictions
from q2_types.sample_data import SampleData
from q2_types.tree import Phylogeny, Rooted

from q2_confluence import __version__
from q2_confluence._examples import (frequency_one_million,
                                     gtdb_tree_collapse_genus,
                                     integrate_at_genus,
                                     merge_metadata_union,
                                     paired_report_example,
                                     profile_reads_example,
                                     summarize_inputs_two_studies)
from q2_confluence._inputs import summarize_inputs
from q2_confluence._integrate import integrate
from q2_confluence._metadata import merge_metadata
from q2_confluence._report import paired_report, summarize_paired_features
from q2_confluence._sylph import (SylphDatabase, SylphDatabaseDirFmt, SylphDbFormat,
                                  SylphTaxonomyFormat, fetch_sylph_database,
                                  profile_reads)
from q2_confluence._table import frequency
from q2_confluence._tree import gtdb_tree_collapse

citations = Citations.load("citations.bib", package="q2_confluence")

plugin = Plugin(
    name="confluence",
    version=__version__,
    website="https://github.com/gregcaporaso/q2-confluence",
    package="q2_confluence",
    description=("Integrate 16S rRNA amplicon and shotgun metagenome taxonomic "
                 "profiles on one GTDB release: collapse each profile to a chosen "
                 "taxonomic level, merge the studies into one feature table with "
                 "provenance, and compare paired samples."),
    short_description="Integrate amplicon and metagenome taxonomic profiles.",
    citations=[citations['Caporaso-Bolyen-2024']]
)

# Provisional actions are registered with deprecated=True and name the plugin
# they are expected to move to, so users see up front what is not stable here.

plugin.methods.register_function(
    function=gtdb_tree_collapse,
    inputs={'tree': Phylogeny[Rooted], 'taxonomy': FeatureData[Taxonomy]},
    parameters={'level': Int % Range(1, None)},
    outputs=[('collapsed_tree', Phylogeny[Rooted])],
    input_descriptions={
        'tree': 'A reference tree whose tips are keyed like the taxonomy, '
                'e.g. the GTDB bac120 tree with genome accessions as tips.',
        'taxonomy': 'The taxonomy of every tip, in the ;-separated rank form.'},
    parameter_descriptions={
        'level': 'The taxonomic level to collapse to, counted as in '
                 'taxa collapse (1 = the first rank; with the seven GTDB ranks '
                 '5 = family, 6 = genus, 7 = species).'},
    output_descriptions={
        'collapsed_tree': 'One tip per taxon at the level, named with the '
                          'taxonomy string collapsed to that level, so the tips '
                          'match the feature IDs of a table collapsed to the '
                          'same level with taxa collapse.'},
    name='Collapse a reference tree to one tip per taxon',
    description=(
        'Group the tips by their taxonomy at the level and replace each group '
        'by its most recent common ancestor, whose branch length becomes its own '
        'length plus the mean distance to the group\'s tips; a group of one tip '
        'is renamed. The action fails if a tip is missing from the taxonomy, a '
        'taxon has no name at the level, or a group is not monophyletic. Run '
        'once per reference release and level; UniFrac tolerates tips absent '
        'from a table, so the same tree serves every table collapsed to that '
        'level. Provisional: expected to move to q2-phylogeny.'),
    deprecated=True,
    examples={'Collapse a GTDB tree to genus': gtdb_tree_collapse_genus},
)

plugin.methods.register_function(
    function=frequency,
    inputs={'table': FeatureTable[RelativeFrequency]},
    parameters={'total': Int % Range(1, None)},
    outputs=[('table', FeatureTable[Frequency])],
    input_descriptions={'table': 'Proportions per sample (each sample sums to 1).'},
    parameter_descriptions={
        'total': 'The count every sample is scaled to. A feature below '
                 '0.5 / total in a sample rounds to zero there, so a total of '
                 '1000000 keeps every proportion of 5e-7 or more.'},
    output_descriptions={
        'table': 'round(proportion * total) per sample; sample totals equal '
                 'the total up to rounding.'},
    name='Rescale proportions to counts',
    description=(
        'Profiles reported as proportions (sylph, MetaPhlAn, TPM) cannot go '
        'through actions that need counts, such as rarefaction, q2-boots, or '
        'the sample classifier. This action multiplies each sample by one '
        'total and rounds, and says so: unlike sapienns frequency, the total '
        'is applied to proportions, not percentages. Provisional: expected to '
        'move to q2-feature-table.'),
    deprecated=True,
    examples={'Scale to one million per sample': frequency_one_million},
)

plugin.visualizers.register_function(
    function=summarize_inputs,
    inputs={'tables': Collection[FeatureTable[Frequency | RelativeFrequency]],
            'metadata': Collection[ImmutableMetadata]},
    parameters={},
    input_descriptions={
        'tables': 'The feature tables to integrate, keyed by a label.',
        'metadata': 'The studies\' sample metadata, keyed by a label '
                    '(import a metadata TSV as ImmutableMetadata).'},
    parameter_descriptions={},
    name='Summarize the sample IDs shared between inputs',
    description=(
        'Before integrating: how many samples each table and metadata input '
        'holds, how many sample IDs each pair shares, and how many IDs fall '
        'into each combination of inputs, with the per-sample membership as a '
        'TSV. Integration pairs samples by identical IDs, so this shows which '
        'samples will pair and which metadata rows have no table.'),
    examples={'Two studies': summarize_inputs_two_studies},
)

plugin.methods.register_function(
    function=merge_metadata,
    inputs={'metadata': Collection[ImmutableMetadata]},
    parameters={'how': Str % Choices(['union', 'intersection'])},
    outputs=[('merged_metadata', ImmutableMetadata)],
    input_descriptions={
        'metadata': 'The studies\' sample metadata, keyed by a label.'},
    parameter_descriptions={
        'how': 'union keeps every sample ID in any input; intersection keeps '
               'only the IDs present in every input (an error if there are '
               'none).'},
    output_descriptions={'merged_metadata': 'One metadata table over the '
                                            'chosen IDs and every column.'},
    name='Merge the metadata of several studies',
    description=(
        'Outer or inner join on sample ID over any number of metadata inputs. '
        'Columns with the same name are combined; a sample with two different '
        'values in such a column is an error. Unlike metadata merge, this takes '
        'a collection and offers the union.'),
    examples={'Union of two studies': merge_metadata_union},
)

plugin.pipelines.register_function(
    function=integrate,
    inputs={'tables': Collection[FeatureTable[Frequency]],
            'taxonomies': Collection[FeatureData[Taxonomy]],
            'tree': Phylogeny[Rooted],
            'metadata': Collection[ImmutableMetadata]},
    parameters={'level': Int % Range(1, None),
                'total': Int % Range(1, None)},
    outputs=[('relative_table', FeatureTable[RelativeFrequency]),
             ('frequency_table', FeatureTable[Frequency]),
             ('taxonomy', FeatureData[Taxonomy]),
             ('integrated_metadata', ImmutableMetadata)],
    input_descriptions={
        'tables': 'One count table per source, keyed by a label (profiles '
                  'reported as proportions go through frequency first).',
        'taxonomies': 'The taxonomy of each table\'s features, same keys, all '
                      'from the same reference release.',
        'tree': 'The reference tree collapsed to the same level with '
                'gtdb-tree-collapse.',
        'metadata': 'The studies\' sample metadata, keyed by a label '
                    '(optional).'},
    parameter_descriptions={
        'level': 'The taxonomic level to collapse to, as in taxa collapse.',
        'total': 'The per-sample total of the frequency table output.'},
    output_descriptions={
        'relative_table': 'Every source collapsed to the level, unnamed features '
                          'dropped, as proportions; sample IDs are <id>-<key>.',
        'frequency_table': 'The same table as counts at the given total.',
        'taxonomy': 'The taxonomy of the merged features (their IDs).',
        'integrated_metadata': 'One row per merged sample: original-sample-id, '
                    'source (the key), unassigned-share (dropped at the level), '
                    'and the union of the studies\' metadata.'},
    name='Integrate profiles at one taxonomic level',
    description=(
        'For each source: taxa collapse to the level, drop the features without '
        'a name at that level, convert to proportions, suffix the sample IDs '
        'with the key. Merge, check every feature is a tip of the tree, rescale '
        'to counts, derive the taxonomy from the feature IDs, and build the '
        'metadata. Samples with the same ID in two sources are a pair, which '
        'paired-report uses.'),
    examples={'Amplicon and metagenome at genus': integrate_at_genus},
)

plugin.register_formats(SylphDbFormat, SylphTaxonomyFormat, SylphDatabaseDirFmt)
plugin.register_semantic_types(SylphDatabase)
plugin.register_artifact_class(
    SylphDatabase, SylphDatabaseDirFmt,
    description='A sylph genome sketch database with the sylph-tax taxonomy '
                'of its genomes.')

plugin.methods.register_function(
    function=fetch_sylph_database,
    inputs={},
    parameters={'release': Str % Choices(['226', '232'])},
    outputs=[('database', SylphDatabase)],
    input_descriptions={},
    parameter_descriptions={'release': 'The GTDB release.'},
    output_descriptions={'database': 'The pre-built database of the GTDB '
                                     'species representatives (c=200, about 14 '
                                     'to 20 GB) with its taxonomy.'},
    name='Fetch a pre-built sylph GTDB database',
    description=(
        'Downloads the sylph authors\' pre-built database of the GTDB species '
        'representatives and the matching sylph-tax taxonomy metadata from '
        'Zenodo. No checksum is published for the database. Provisional: '
        'expected to move to MOSHPIT.'),
    deprecated=True,
)

plugin.methods.register_function(
    function=profile_reads,
    inputs={'reads': SampleData[PairedEndSequencesWithQuality],
            'database': SylphDatabase},
    parameters={'threads': Int % Range(1, None), 'estimate_unknown': Bool},
    outputs=[('taxonomic_abundance', FeatureTable[RelativeFrequency]),
             ('sequence_abundance', FeatureTable[RelativeFrequency]),
             ('taxonomy', FeatureData[Taxonomy]),
             ('stats', ImmutableMetadata)],
    input_descriptions={
        'reads': 'Paired-end reads per sample (quality-filtered, host reads '
                 'removed).',
        'database': 'A sylph database with its taxonomy.'},
    parameter_descriptions={
        'threads': 'Threads for sylph sketch and profile.',
        'estimate_unknown': 'Estimate the share of reads from no genome in the '
                            'database (sylph profile -u); it is reported per '
                            'sample as unknown-share and the sequence abundance '
                            'is renormalized over the detected genomes.'},
    output_descriptions={
        'taxonomic_abundance': 'Per sample, the share of each detected genome '
                               'corrected for genome size (roughly a share of '
                               'cells; the measure to integrate).',
        'sequence_abundance': 'Per sample, the share of reads from each '
                              'detected genome.',
        'taxonomy': 'The taxonomy of every detected genome, keyed by accession, '
                    'seven ranks joined with ;.',
        'stats': 'Per sample: genomes and genera detected, median ANI, '
                 'unknown share.'},
    name='Profile reads against a sylph database',
    description=(
        'sylph sketch on every sample, sylph profile against the database, and '
        'sylph-tax taxprof with the database\'s taxonomy, parsed into '
        'genome-keyed tables. A sample in which sylph detects nothing keeps an '
        'all-zero column. Provisional: expected to move to MOSHPIT.'),
    deprecated=True,
    examples={'Profile three samples against a two-genome database':
              profile_reads_example},
)

plugin.visualizers.register_function(
    function=summarize_paired_features,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency],
            'importance': FeatureData[Importance],
            'predictions': SampleData[ClassifierPredictions]},
    parameters={'metadata': Metadata},
    input_descriptions={
        'table': 'The integrated table of the paired samples.',
        'importance': 'Feature importances from the random forest.',
        'predictions': 'Cross-validated predictions of the source.'},
    parameter_descriptions={
        'metadata': 'The integrated metadata (columns source and '
                    'original-sample-id).'},
    name='Summarize the features that separate the sources',
    description=(
        'Per feature: random-forest importance and rank, each source\'s mean '
        'relative abundance and prevalence, the log2 ratio between the first '
        'two sources, and whether the feature is seen by one source only; '
        'importance shares by direction and the cross-validated accuracy. '
        'The full table is downloadable as a TSV.'),
)

plugin.pipelines.register_function(
    function=paired_report,
    inputs={'table': FeatureTable[Frequency],
            'tree': Phylogeny[Rooted],
            'metadata': ImmutableMetadata},
    parameters={'n_boots': Int % Range(1, None),
                'permutations': Int % Range(1, None),
                'cv': Int % Range(2, None),
                'n_estimators': Int % Range(1, None),
                'random_state': Int,
                'n_jobs': Int % Range(1, None)},
    outputs=[('report', Visualization),
             ('distance_matrices', Collection[DistanceMatrix]),
             ('feature_importance', FeatureData[Importance])],
    input_descriptions={
        'table': 'The frequency table from integrate.',
        'tree': 'The reference tree collapsed to the same level.',
        'metadata': 'The metadata from integrate (source and '
                    'original-sample-id).'},
    parameter_descriptions={
        'n_boots': 'Bootstrap iterations for q2-boots core-metrics.',
        'permutations': 'PERMANOVA permutations.',
        'cv': 'Folds of the nested cross-validation.',
        'n_estimators': 'Trees in the random forest.',
        'random_state': 'Seed for the bootstraps and the forest.',
        'n_jobs': 'Threads for the random forest.'},
    output_descriptions={
        'report': 'One report: the bootstrapped ordination, PERMANOVA by source '
                  'for every distance metric, the random forest\'s confusion '
                  'matrix and heatmap, and the per-feature summary.',
        'distance_matrices': 'The bootstrapped distance matrices.',
        'feature_importance': 'The random forest\'s feature importances.'},
    name='Compare the sources on the samples every source profiled',
    description=(
        'Keeps the samples whose original ID is present in every source, runs '
        'q2-boots core-metrics with the tree at the smallest sample total, '
        'PERMANOVA by source on every distance matrix, a random forest '
        'classifying source with nested cross-validation, and the per-feature '
        'summary, collated into one report.'),
    examples={'Two sources on six paired samples': paired_report_example},
)
