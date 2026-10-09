# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Plugin
from q2_types.feature_table import FeatureTable, Frequency
from la_confluence import __version__
from la_confluence._methods import duplicate_table

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
    # The plugin-level citation of 'Caporaso-Bolyen-2024' is provided as
    # an example. You can replace this with citations to other references
    # in citations.bib.
    citations=[citations['Caporaso-Bolyen-2024']]
)

plugin.methods.register_function(
    function=duplicate_table,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[('new_table', FeatureTable[Frequency])],
    input_descriptions={'table': 'The feature table to be duplicated.'},
    parameter_descriptions={},
    output_descriptions={'new_table': 'The duplicated feature table.'},
    name='Duplicate table',
    description=("Create a copy of a feature table with a new uuid. "
                 "This is for demonstration purposes only. "),
    citations=[]
)
