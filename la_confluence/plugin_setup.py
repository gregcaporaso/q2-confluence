# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Plugin

from la_confluence import __version__

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
