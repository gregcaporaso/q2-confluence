# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Collection, Int, Metadata, Plugin, Range,
                           Visualization)
from q2_types.feature_data import FeatureData, Importance, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.metadata import ImmutableMetadata
from q2_types.tree import Phylogeny, Rooted

from q2_confluence import __version__
from q2_confluence._examples import (_relative_frequency_example,
                                     _gtdb_tree_collapse_genus,
                                     _integrate_genomes,
                                     _paired_report_example,
                                     _summarize_inputs_two_studies)
from q2_confluence._formats import GTDBNewickFormat, _gtdb_newick_to_newick
from q2_confluence._inputs import summarize_inputs
from q2_confluence._integrate import integrate
from q2_confluence._report import paired_report, _summarize_paired_features
from q2_confluence._table import relative_frequency
from q2_confluence._tree import GTDB_RANKS, gtdb_tree_collapse

plugin = Plugin(
    name="confluence",
    version=__version__,
    website="https://github.com/gregcaporaso/q2-confluence",
    package="q2_confluence",
    description=("Integrate 16S rRNA amplicon and metagenome taxonomic "
                 "profiles in the feature space of one GTDB release's genomes: "
                 "filter each profile to the reference tree, merge the studies "
                 "into one feature table with provenance, and compare paired "
                 "samples."),
    short_description="Integrate amplicon and metagenome taxonomic profiles.",
)

# Provisional actions name the plugin they are expected to move to in their
# description; the README's alpha note covers the plugin as a whole.

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
        'level. Provisional: expected to move to q2-phylogeny.'),
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
        'feature-table relative-frequency.'),
    examples={'Renormalize after dropping a feature': _relative_frequency_example},
)

plugin.visualizers.register_function(
    function=summarize_inputs,
    inputs={'tables': Collection[FeatureTable[Frequency | RelativeFrequency]],
            'tree': Phylogeny[Rooted],
            'metadata': Collection[ImmutableMetadata]},
    parameters={},
    input_descriptions={
        'tables': 'The tables integrate will take, keyed by study.',
        'tree': 'The reference tree integrate will take; a study\'s features '
                'that are not its tips count as dropped.',
        'metadata': 'The studies\' sample metadata, keyed like the tables '
                    '(optional, any subset of the keys).'},
    parameter_descriptions={},
    name='Summarize the inputs of integrate',
    description=(
        'Before integrating: per study, how many samples its table and its '
        'metadata hold and share, and how many of its features are tips of '
        'the tree (what integrate keeps); across studies, the sample IDs the '
        'tables share (integrate pairs samples by identical IDs) and the '
        'feature IDs the tables and the tree share; the per-sample '
        'membership as a TSV.'
        ' EXPERIMENTAL: this is AI-generated code; human review pending.'),
    examples={'Two studies': _summarize_inputs_two_studies},
)

plugin.pipelines.register_function(
    function=integrate,
    inputs={'tables': Collection[FeatureTable[Frequency | RelativeFrequency]],
            'tree': Phylogeny[Rooted],
            'metadata': Collection[ImmutableMetadata]},
    parameters={},
    outputs=[('integrated_table', FeatureTable[RelativeFrequency]),
             ('integrated_metadata', ImmutableMetadata)],
    input_descriptions={
        'tables': 'One table per source, counts or proportions, keyed by a '
                  'label; every table is keyed like the tree\'s tips (GTDB '
                  'genome accessions, from closed-reference clustering on '
                  'amplicon data and from sylph profile-to-tables on '
                  'metagenome data). One table alone works too, to prepare '
                  'for manual integration.',
        'tree': 'The reference tree as imported (the GTDB tree through '
                'GTDBNewickFormat); features that are not its tips are '
                'dropped from every table.',
        'metadata': 'The studies\' sample metadata, keyed like the tables '
                    '(optional, any subset of the keys; a study file shared '
                    'by two sources is given under each key). Each study\'s '
                    'IDs are suffixed with its key, the studies are merged '
                    'with metadata merge and joined to the output metadata.'},
    parameter_descriptions={},
    output_descriptions={
        'integrated_table': 'Every source as proportions over the features '
                            'kept, merged; sample IDs are <id>-<key>.',
        'integrated_metadata': 'One row per merged sample: '
                               'original-sample-id, study-id (the key), '
                               'the source\'s features-before and '
                               'features-after the tree filter and their '
                               'quotient features-kept-share, and the columns '
                               'of the sample\'s study metadata.'},
    name='Integrate profiles in the feature space of the reference genomes',
    description=(
        'For each source: drop the features that are not tips of the tree '
        '(phylogeny filter-table; counted in the metadata, with a warning '
        'when more than a tenth goes), renormalize to proportions and suffix '
        'the sample IDs with the key. Merge the sources and build the '
        'metadata. The taxonomy of the merged features is the reference\'s '
        'own (GTDB\'s taxonomy file imported as FeatureData[Taxonomy]), which '
        'any downstream action takes directly. Samples with the same ID in '
        'two sources are considered a pair in paired-report.'),
    examples={'Amplicon and metagenome in the genome space':
              _integrate_genomes},
)

plugin.methods.register_function(
    function=_summarize_paired_features,
    inputs={'table': FeatureTable[RelativeFrequency],
            'importance': FeatureData[Importance]},
    parameters={'metadata': Metadata},
    outputs=[('feature_summary', ImmutableMetadata)],
    input_descriptions={
        'table': 'The integrated table of the paired samples, from exactly '
                 'two studies.',
        'importance': 'Feature importances from the random forest.'},
    parameter_descriptions={
        'metadata': 'The integrated metadata (columns study-id and '
                    'original-sample-id).'},
    output_descriptions={
        'feature_summary': 'One row per feature: importance and rank, '
                           'mean-<study-id> and prevalence-<study-id> for '
                           'each of the two studies, log2-<second>-over-<first>, '
                           'and direction (shared or <study-id>-only).'},
    name='Summarize the features that separate the sources',
    description=(
        'Per feature: random-forest importance and rank, each study\'s mean '
        'relative abundance and prevalence, the log2 ratio of the second '
        'study over the first, and whether the feature is seen by one study '
        'only; the table must hold exactly two studies. '
        'Metadata rather than a visualization, so it can be tabulated, '
        'filtered with where clauses, or joined to a taxonomy.'),
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
        'table': 'The integrated table from integrate; its samples must belong '
                 'to exactly two studies (feature-table filter-samples '
                 'keeps two).',
        'tree': 'The reference tree integrate was run with.',
        'metadata': 'The metadata from integrate (study-id and '
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
        'report': 'One report: per metric, PERMANOVA by study-id and a PCoA '
                  'scatter plot; the random forest\'s confusion matrix; and '
                  'the per-feature summary tabulated (every feature, '
                  'searchable and sortable, with the taxonomy when given).',
        'distance_matrices': 'The distance matrices, keyed braycurtis and '
                             'weighted_unifrac.',
        'feature_importance': 'The random forest\'s feature importances.',
        'feature_summary': 'The per-feature summary as metadata.'},
    name='Compare the sources on the samples every source profiled',
    description=(
        'Keeps the samples whose original ID is present in both studies and '
        'computes Bray-Curtis and weighted UniFrac on the proportions as they '
        'are (no rarefaction, so no presence-based metric), with PERMANOVA by '
        'study-id and a PCoA scatter plot for each; then a random forest classifying '
        'study-id with nested cross-validation and the per-feature summary, '
        'collated into one report.'),
    examples={'Two sources on six paired samples': _paired_report_example},
)
