# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import numpy as np
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_confluence import _examples as ex


class IntegrateTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def setUp(self):
        # TODO: this is AI-generated code; human review pending
        super().setUp()
        self.integrate = self.plugin.pipelines['integrate']
        self.tree = ex._reference_tree_factory()
        taxonomy = ex._reference_taxonomy_factory()
        self.tables = {'amplicon': ex._amplicon_genome_table_factory(),
                       'metagenome': ex._metagenome_genome_table_factory()}
        self.taxonomies = {'amplicon': taxonomy, 'metagenome': taxonomy}
        self.metadata = {'amplicon': ex._amplicon_metadata_factory(),
                         'metagenome': ex._metagenome_metadata_factory()}

    def test_genomes(self):
        # TODO: this is AI-generated code; human review pending
        rel, taxonomy, md = self.integrate(self.tables, self.taxonomies,
                                           self.tree, metadata=self.metadata)
        rel = rel.view(biom.Table)
        self.assertEqual(sorted(rel.ids(axis='observation')),
                         ['GB_GCA_3', 'RS_GCF_1', 'RS_GCF_2', 'RS_GCF_4'])
        self.assertEqual(sorted(rel.ids()),
                         ['sample-1-amplicon', 'sample-2-amplicon',
                          'sample-2-metagenome', 'sample-3-amplicon',
                          'sample-3-metagenome'])
        # nothing is collapsed or dropped: sample-1 amplicon = 12 + 3 + 1 counts
        self.assertAlmostEqual(rel.get_value_by_ids('RS_GCF_1', 'sample-1-amplicon'),
                               12 / 16)
        # the metagenome source came in as proportions
        self.assertAlmostEqual(rel.get_value_by_ids('RS_GCF_2', 'sample-2-metagenome'),
                               0.4)
        tax = taxonomy.view(pd.Series)
        self.assertEqual(sorted(tax.index), sorted(rel.ids(axis='observation')))
        self.assertEqual(tax['GB_GCA_3'], ex._TAXONOMY['GB_GCA_3'])
        md = md.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(md.loc['sample-1-amplicon', 'original-sample-id'], 'sample-1')
        self.assertEqual(md.loc['sample-1-amplicon', 'profile-type'], 'amplicon')
        self.assertEqual(md.loc['sample-1-amplicon', 'site'], 'compost')
        self.assertAlmostEqual(md.loc['sample-3-metagenome', 'read-pairs-millions'],
                               38.7)

    def test_without_metadata(self):
        # TODO: this is AI-generated code; human review pending
        _, _, md = self.integrate(self.tables, self.taxonomies, self.tree)
        cols = md.view(qiime2.Metadata).to_dataframe().columns
        self.assertEqual(sorted(cols), ['original-sample-id', 'profile-type'])

    def test_missing_taxonomy(self):
        # TODO: this is AI-generated code; human review pending
        partial = ex._taxonomy_factory({'RS_GCF_1': ex._TAXONOMY['RS_GCF_1']})
        with self.assertRaisesRegex(ValueError, 'features have no taxonomy'):
            self.integrate(self.tables,
                           {'amplicon': partial, 'metagenome': partial},
                           self.tree)

    def test_not_tips(self):
        # TODO: this is AI-generated code; human review pending
        # a table keyed by unprefixed accessions is not in the tree's space
        table = biom.Table(np.array([[0.6, 0.1], [0.4, 0.9]]),
                           ['GCF_1', 'GCF_2'], ['sample-2', 'sample-3'])
        unprefixed = qiime2.Artifact.import_data(
            'FeatureTable[RelativeFrequency]', table)
        taxonomy = ex._taxonomy_factory({'GCF_1': ex._TAXONOMY['RS_GCF_1'],
                                         'GCF_2': ex._TAXONOMY['RS_GCF_2']})
        with self.assertRaisesRegex(ValueError, 'not tips of the tree'):
            self.integrate({'amplicon': self.tables['amplicon'],
                            'metagenome': unprefixed},
                           {'amplicon': self.taxonomies['amplicon'],
                            'metagenome': taxonomy},
                           self.tree)

    def test_key_mismatch(self):
        # TODO: this is AI-generated code; human review pending
        with self.assertRaisesRegex(ValueError, 'same keys'):
            self.integrate(self.tables, {'amplicon': self.taxonomies['amplicon']},
                           self.tree)
