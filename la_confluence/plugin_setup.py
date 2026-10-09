# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Int, Plugin, Range
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.tree import Phylogeny, Rooted

from la_confluence import __version__
from la_confluence._examples import (frequency_one_million,
                                     gtdb_tree_collapse_genus)
from la_confluence._table import frequency
from la_confluence._tree import gtdb_tree_collapse

citations = Citations.load("citations.bib", package="la_confluence")

plugin = Plugin(
    name="confluence",
    version=__version__,
    website="https://github.com/gregcaporaso/la-confluence",
    package="la_confluence",
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
