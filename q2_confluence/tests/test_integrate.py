# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_confluence import _examples as ex

BLAUTIA = ('d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
           'f__Lachnospiraceae;g__Blautia')
BACILLUS = ('d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;'
            'g__Bacillus')


class IntegrateTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def setUp(self):
        super().setUp()
        self.integrate = self.plugin.pipelines['integrate']
        collapse = self.plugin.methods['gtdb_tree_collapse']
        self.genus_tree, = collapse(ex.reference_tree_factory(),
                                    ex.reference_taxonomy_factory(), 6)
        self.tables = {'amplicon': ex.amplicon_table_three_asvs_factory(),
                       'metagenome': ex.metagenome_counts_factory()}
        self.taxonomies = {'amplicon': ex.amplicon_taxonomy_factory(),
                           'metagenome': ex.metagenome_taxonomy_factory()}
        self.metadata = {'amplicon': ex.amplicon_metadata_factory(),
                         'metagenome': ex.metagenome_metadata_factory()}

    def test_genus(self):
        rel, counts, taxonomy, md = self.integrate(
            self.tables, self.taxonomies, self.genus_tree, 6, 1000,
            metadata=self.metadata)
        rel = rel.view(biom.Table)
        self.assertEqual(sorted(rel.ids(axis='observation')), [BACILLUS, BLAUTIA])
        self.assertEqual(sorted(rel.ids()),
                         ['sample-1-amplicon', 'sample-2-amplicon',
                          'sample-2-metagenome', 'sample-3-amplicon',
                          'sample-3-metagenome'])
        # asv-3 (unnamed at genus) is dropped and the sample renormalized:
        # sample-1 amplicon = 12 Blautia, 3 Bacillus
        self.assertAlmostEqual(rel.get_value_by_ids(BLAUTIA, 'sample-1-amplicon'),
                               12 / 15)
        counts = counts.view(biom.Table)
        self.assertEqual(counts.get_value_by_ids(BLAUTIA, 'sample-1-amplicon'), 800)
        self.assertEqual(counts.get_value_by_ids(BLAUTIA, 'sample-2-metagenome'), 600)
        tax = taxonomy.view(pd.Series)
        self.assertEqual(tax[BLAUTIA], BLAUTIA)
        md = md.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(md.loc['sample-1-amplicon', 'original-sample-id'], 'sample-1')
        self.assertEqual(md.loc['sample-1-amplicon', 'source'], 'amplicon')
        self.assertAlmostEqual(md.loc['sample-1-amplicon', 'unassigned-share'],
                               1 / 16)
        self.assertEqual(md.loc['sample-1-amplicon', 'site'], 'compost')
        self.assertAlmostEqual(md.loc['sample-3-metagenome', 'read-pairs-millions'],
                               38.7)
        self.assertEqual(md.loc['sample-2-metagenome', 'unassigned-share'], 0)

    def test_without_metadata(self):
        _, _, _, md = self.integrate(self.tables, self.taxonomies,
                                     self.genus_tree, 6, 1000)
        cols = md.view(qiime2.Metadata).to_dataframe().columns
        self.assertEqual(sorted(cols),
                         ['original-sample-id', 'source', 'unassigned-share'])

    def test_tree_at_another_level(self):
        family_tree, = self.plugin.methods['gtdb_tree_collapse'](
            ex.reference_tree_factory(), ex.reference_taxonomy_factory(), 5)
        with self.assertRaisesRegex(ValueError, 'not tips of the tree'):
            self.integrate(self.tables, self.taxonomies, family_tree, 6, 1000)

    def test_key_mismatch(self):
        with self.assertRaisesRegex(ValueError, 'same keys'):
            self.integrate(self.tables, {'amplicon': self.taxonomies['amplicon']},
                           self.genus_tree, 6, 1000)
