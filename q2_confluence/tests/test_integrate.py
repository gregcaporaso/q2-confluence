# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import numpy as np
import qiime2
from qiime2.plugin.testing import TestPluginBase
from rachis.core.exceptions import RachisWarning

from q2_confluence import _examples as ex


class IntegrateTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def setUp(self):
        # TODO: this is AI-generated code; human review pending
        super().setUp()
        self.integrate = self.plugin.pipelines['integrate']
        self.tree = ex._reference_tree_factory()
        self.tables = {'amplicon': ex._amplicon_genome_table_factory(),
                       'metagenome': ex._metagenome_genome_table_factory()}
        self.metadata = {'amplicon': ex._amplicon_metadata_factory(),
                         'metagenome': ex._metagenome_metadata_factory()}

    def test_genomes(self):
        # TODO: this is AI-generated code; human review pending
        rel, md = self.integrate(self.tables, self.tree, metadata=self.metadata)
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
        md = md.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(list(md.columns)[:5],
                         ['original-sample-id', 'study-id', 'features-before',
                          'features-after', 'features-kept-share'])
        self.assertEqual(md.loc['sample-1-amplicon', 'original-sample-id'], 'sample-1')
        self.assertEqual(md.loc['sample-1-amplicon', 'study-id'], 'amplicon')
        self.assertEqual(md.loc['sample-1-amplicon', 'features-before'], 3)
        self.assertEqual(md.loc['sample-1-amplicon', 'features-after'], 3)
        self.assertEqual(md.loc['sample-1-amplicon', 'features-kept-share'], 1.0)
        # each merged sample carries its own study's columns: the studies
        # share IDs and the column site, which is fine once the IDs are
        # suffixed per source
        self.assertEqual(md.loc['sample-1-amplicon', 'site'], 'compost')
        self.assertEqual(md.loc['sample-3-amplicon', 'site'], 'soil')
        self.assertEqual(md.loc['sample-3-metagenome', 'site'], 'soil')
        self.assertAlmostEqual(md.loc['sample-3-metagenome', 'read-pairs-millions'],
                               38.7)
        self.assertTrue(np.isnan(md.loc['sample-3-amplicon', 'read-pairs-millions']))

    def test_metadata_key_not_a_table_key(self):
        # TODO: this is AI-generated code; human review pending
        with self.assertRaisesRegex(ValueError, r"metadata keys \['study'\] are "
                                    'not table keys'):
            self.integrate(self.tables, self.tree,
                           metadata={'study': self.metadata['amplicon']})

    def test_without_metadata(self):
        # TODO: this is AI-generated code; human review pending
        _, md = self.integrate(self.tables, self.tree)
        cols = md.view(qiime2.Metadata).to_dataframe().columns
        self.assertEqual(sorted(cols), ['features-after', 'features-before',
                                        'features-kept-share',
                                        'original-sample-id', 'study-id'])

    def test_single_table(self):
        # TODO: this is AI-generated code; human review pending
        rel, md = self.integrate({'amplicon': self.tables['amplicon']}, self.tree,
                                 metadata={'amplicon': self.metadata['amplicon']})
        rel = rel.view(biom.Table)
        self.assertEqual(sorted(rel.ids()), ['sample-1-amplicon',
                                             'sample-2-amplicon',
                                             'sample-3-amplicon'])
        self.assertTrue(np.allclose(rel.sum(axis='sample'), 1.0))
        md = md.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(md.loc['sample-2-amplicon', 'site'], 'compost')

    def test_features_not_in_tree_dropped_and_counted(self):
        # TODO: this is AI-generated code; human review pending
        # two of three features are not tips (an archaeon against the
        # bacterial tree, say): dropped, renormalized, counted, and warned
        # about since more than a tenth went
        table = biom.Table(np.array([[0.5, 0.2], [0.3, 0.5], [0.2, 0.3]]),
                           ['RS_GCF_1', 'GB_GCA_archaeon', 'GB_GCA_other'],
                           ['sample-2', 'sample-3'])
        extra = qiime2.Artifact.import_data('FeatureTable[RelativeFrequency]',
                                            table)
        with self.assertWarnsRegex(RachisWarning, 'metagenome: 2 of its 3 '
                                   'features .67%. are not tips'):
            rel, md = self.integrate({'amplicon': self.tables['amplicon'],
                                      'metagenome': extra}, self.tree)
        rel = rel.view(biom.Table)
        self.assertEqual(sorted(rel.ids(axis='observation')),
                         ['GB_GCA_3', 'RS_GCF_1', 'RS_GCF_4'])
        self.assertAlmostEqual(rel.get_value_by_ids('RS_GCF_1', 'sample-2-metagenome'),
                               1.0)
        md = md.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(md.loc['sample-2-metagenome', 'features-before'], 3)
        self.assertEqual(md.loc['sample-2-metagenome', 'features-after'], 1)
        self.assertAlmostEqual(md.loc['sample-2-metagenome', 'features-kept-share'],
                               1 / 3)
        self.assertEqual(md.loc['sample-2-amplicon', 'features-kept-share'], 1.0)

    def test_no_feature_in_tree(self):
        # TODO: this is AI-generated code; human review pending
        # a table keyed by unprefixed accessions is not in the tree's space
        table = biom.Table(np.array([[0.6, 0.1], [0.4, 0.9]]),
                           ['GCF_1', 'GCF_2'], ['sample-2', 'sample-3'])
        unprefixed = qiime2.Artifact.import_data(
            'FeatureTable[RelativeFrequency]', table)
        with self.assertRaisesRegex(ValueError, 'metagenome: none of its 2 '
                                    'features is a tip of the tree'):
            self.integrate({'amplicon': self.tables['amplicon'],
                            'metagenome': unprefixed}, self.tree)
