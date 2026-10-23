# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Choices, Collection, Int, Metadata, Plugin, Range,
                           Str, Visualization)
from q2_types.feature_data import FeatureData, Importance, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.metadata import ImmutableMetadata
from q2_types.tree import Phylogeny, Rooted

from q2_confluence import __version__
from q2_confluence._examples import (_relative_frequency_example,
                                     _gtdb_tree_collapse_genus,
                                     _integrate_genomes,
                                     _merge_metadata_union,
                                     _paired_report_example,
                                     _summarize_paired_features_example,
                                     _summarize_inputs_two_studies)
from q2_confluence._formats import GTDBNewickFormat, _gtdb_newick_to_newick
from q2_confluence._inputs import summarize_inputs
from q2_confluence._integrate import integrate
from q2_confluence._metadata import merge_metadata
from q2_confluence._report import paired_report, summarize_paired_features
from q2_confluence._table import relative_frequency
from q2_confluence._tree import GTDB_RANKS, gtdb_tree_collapse

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
)

# Provisional actions are registered with deprecated=True and name the plugin
# they are expected to move to, so users see up front what is not stable here.

# GTDB's released trees import as Phylogeny[Rooted] through this format
# (qiime tools import --input-format GTDBNewickFormat); the transformer drops
# the root label UniFrac cannot parse. Provisional: expected to move to
# q2-types or q2-phylogeny.
plugin.register_formats(GTDBNewickFormat)
plugin.register_transformer(_gtdb_newick_to_newick)

plugin.methods.register_function(
    function=gtdb_tree_collapse,
    inputs={'tree': Phylogeny[Rooted], 'taxonomy': FeatureData[Taxonomy]},
    parameters={'level': Int % Range(1, GTDB_RANKS, inclusive_end=True)},
    outputs=[('collapsed_tree', Phylogeny[Rooted]),
             ('collapse_log', ImmutableMetadata)],
    input_descriptions={
        'tree': 'A reference tree whose tips are keyed like the taxonomy: '
                'the GTDB tree imported through GTDBNewickFormat, which quotes '
                'the genome accessions so that their underscores are kept.',
        'taxonomy': 'The taxonomy of every tip, in the ;-separated rank form.'},
    parameter_descriptions={
        'level': 'The taxonomic level to collapse to, counted as in '
                 'taxa collapse (1 = the first rank; with the seven GTDB ranks '
                 '5 = family, 6 = genus, 7 = species).'},
    output_descriptions={
        'collapsed_tree': 'One tip per taxon at the level, named with the '
                          'taxonomy string collapsed to that level, so the tips '
                          'match the feature IDs of a table collapsed to the '
                          'same level with taxa collapse.',
        'collapse_log': 'One row per tip of the collapsed tree, keyed by the '
                        'tip name (so it serves as feature metadata for a '
                        'table collapsed to the same level): original-name, '
                        'the name the node had in the input tree; genomes, '
                        'the number of genomes it stands for; genome-ids, '
                        'their IDs, comma-separated.'},
    name='Collapse a reference tree to one tip per taxon',
    description=(
        'Group the tips by their taxonomy at the level and shear each group off '
        'at its most recent common ancestor, which becomes the tip and keeps '
        'its own branch length; a group of one tip is renamed; the other '
        'internal node labels are kept. The action fails if a tip is missing '
        'from the taxonomy, a '
        'taxon has no name at the level, or a group is not monophyletic. Run '
        'once per reference release and level; UniFrac tolerates tips absent '
        'from a table, so the same tree serves every table collapsed to that '
        'level. Provisional: expected to move to q2-phylogeny.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    deprecated=True,
    examples={'Collapse a GTDB tree to genus': _gtdb_tree_collapse_genus},
)

plugin.methods.register_function(
    function=relative_frequency,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency]},
    parameters={},
    outputs=[('relative_frequency_table', FeatureTable[RelativeFrequency])],
    input_descriptions={
        'table': 'Counts, or proportions whose samples no longer sum to 1 '
                 '(after a feature filter).'},
    parameter_descriptions={},
    output_descriptions={
        'relative_frequency_table': 'Each sample divided by its total; an '
                                    'empty sample stays empty.'},
    name='Proportions from counts or from proportions',
    description=(
        'feature-table relative-frequency for a table that may already hold '
        'proportions: each sample is divided by its total, so a relative '
        'table from which features were dropped sums to 1 again. integrate '
        'uses it on every source. Provisional: expected to be folded into '
        'feature-table relative-frequency.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    deprecated=True,
    examples={'Renormalize after dropping a feature': _relative_frequency_example},
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
        'samples will pair and which metadata rows have no table.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    examples={'Two studies': _summarize_inputs_two_studies},
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
        'a collection and offers the union.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    examples={'Union of two studies': _merge_metadata_union},
)

plugin.pipelines.register_function(
    function=integrate,
    inputs={'tables': Collection[FeatureTable[Frequency | RelativeFrequency]],
            'taxonomies': Collection[FeatureData[Taxonomy]],
            'tree': Phylogeny[Rooted],
            'metadata': Collection[ImmutableMetadata]},
    parameters={},
    outputs=[('relative_table', FeatureTable[RelativeFrequency]),
             ('taxonomy', FeatureData[Taxonomy]),
             ('integrated_metadata', ImmutableMetadata)],
    input_descriptions={
        'tables': 'One table per source, counts or proportions, keyed by a '
                  'label; every table is keyed by the tree\'s tips (GTDB '
                  'genome accessions, from closed-reference clustering on the '
                  'amplicon side and from sylph profile-to-tables on the '
                  'metagenome side).',
        'taxonomies': 'The taxonomy of each table\'s features, same keys, all '
                      'from the same reference release.',
        'tree': 'The reference tree as imported (the GTDB tree through '
                'GTDBNewickFormat).',
        'metadata': 'The studies\' sample metadata, keyed by a label '
                    '(optional).'},
    parameter_descriptions={},
    output_descriptions={
        'relative_table': 'Every source as proportions, merged; sample IDs '
                          'are <id>-<key>.',
        'taxonomy': 'The input taxonomies restricted to the merged features '
                    '(the first source\'s string when two differ).',
        'integrated_metadata': 'One row per merged sample: original-sample-id, '
                    'profile-type (the key), and the union of the studies\' '
                    'metadata.'},
    name='Integrate profiles in the feature space of the reference genomes',
    description=(
        'For each source: renormalize to proportions and suffix the sample '
        'IDs with the key. Merge, check every feature is a tip of the tree, '
        'take the taxonomy from the inputs, and build the metadata. Samples '
        'with the same ID in two sources are a pair, which paired-report '
        'uses.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    examples={'Amplicon and metagenome in the genome space':
              _integrate_genomes},
)

plugin.methods.register_function(
    function=summarize_paired_features,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency],
            'importance': FeatureData[Importance]},
    parameters={'metadata': Metadata},
    outputs=[('feature_summary', ImmutableMetadata)],
    input_descriptions={
        'table': 'The integrated table of the paired samples.',
        'importance': 'Feature importances from the random forest.'},
    parameter_descriptions={
        'metadata': 'The integrated metadata (columns profile-type and '
                    'original-sample-id).'},
    output_descriptions={
        'feature_summary': 'One row per feature: importance and rank, '
                           'mean-<profile-type> and prevalence-<profile-type> '
                           'for every profile type, log2-<second>-over-<first>, '
                           'and direction (shared or <profile-type>-only).'},
    name='Summarize the features that separate the sources',
    description=(
        'Per feature: random-forest importance and rank, each source\'s mean '
        'relative abundance and prevalence, the log2 ratio between the first '
        'two sources, and whether the feature is seen by one source only. '
        'Metadata rather than a visualization, so it can be tabulated, '
        'filtered with where clauses, or joined to a taxonomy.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    examples={'Three features of two sources': _summarize_paired_features_example},
)

plugin.pipelines.register_function(
    function=paired_report,
    inputs={'table': FeatureTable[RelativeFrequency],
            'tree': Phylogeny[Rooted],
            'metadata': ImmutableMetadata,
            'taxonomy': FeatureData[Taxonomy]},
    parameters={'permutations': Int % Range(1, None),
                'cv': Int % Range(2, None),
                'n_estimators': Int % Range(1, None),
                'random_state': Int,
                'n_jobs': Int % Range(1, None)},
    outputs=[('report', Visualization),
             ('distance_matrices', Collection[DistanceMatrix]),
             ('feature_importance', FeatureData[Importance]),
             ('feature_summary', ImmutableMetadata)],
    input_descriptions={
        'table': 'The relative table from integrate.',
        'tree': 'The reference tree integrate was run with.',
        'metadata': 'The metadata from integrate (profile-type and '
                    'original-sample-id).',
        'taxonomy': 'The taxonomy from integrate (optional): tabulated beside '
                    'the per-feature summary, so a genome ID shows its '
                    'lineage.'},
    parameter_descriptions={
        'permutations': 'PERMANOVA permutations.',
        'cv': 'Folds of the nested cross-validation.',
        'n_estimators': 'Trees in the random forest.',
        'random_state': 'Seed for the forest.',
        'n_jobs': 'Threads for the random forest.'},
    output_descriptions={
        'report': 'One report: per metric, PERMANOVA by profile type and a PCoA '
                  'scatter plot; the random forest\'s confusion matrix; and '
                  'the per-feature summary tabulated (every feature, '
                  'searchable and sortable, with the taxonomy when given).',
        'distance_matrices': 'The distance matrices, keyed by metric.',
        'feature_importance': 'The random forest\'s feature importances.',
        'feature_summary': 'The per-feature summary (summarize-paired-features) '
                           'as metadata.'},
    name='Compare the sources on the samples every source profiled',
    description=(
        'Keeps the samples whose original ID is present in every source and '
        'computes Bray-Curtis, Jaccard, unweighted and weighted UniFrac on the '
        'proportions as they are (no rarefaction), with PERMANOVA by profile type '
        'and a PCoA scatter plot for each; then a random forest classifying '
        'profile type with nested cross-validation and the per-feature summary, '
        'collated into one report.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    examples={'Two sources on six paired samples': _paired_report_example},
)
