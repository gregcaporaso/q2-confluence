# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_confluence import _examples as ex


class PairedReportTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def setUp(self):
        super().setUp()
        collapse = self.plugin.methods['gtdb_tree_collapse']
        self.tree, = collapse(ex.reference_tree_factory(),
                              ex.reference_taxonomy_factory(), 6)
        _, self.table, _, self.metadata = self.plugin.pipelines['integrate'](
            {'amplicon': ex.paired_amplicon_table_factory(),
             'metagenome': ex.paired_metagenome_table_factory()},
            {'amplicon': ex.paired_amplicon_taxonomy_factory(),
             'metagenome': ex.paired_metagenome_taxonomy_factory()},
            self.tree, 6, 1000, metadata={'study': ex.paired_metadata_factory()})

    def test_report(self):
        report, dms, importance = self.plugin.pipelines['paired_report'](
            self.table, self.tree, self.metadata, n_boots=2, permutations=99,
            cv=3, n_estimators=50, random_state=1)
        self.assertEqual(sorted(dms.keys()),
                         ['braycurtis', 'jaccard', 'unweighted_unifrac',
                          'weighted_unifrac'])
        imp = importance.view(pd.DataFrame)
        self.assertEqual(len(imp), 3)
        with tempfile.TemporaryDirectory() as d:
            report.export_data(d)
            names = {f for _, _, files in os.walk(d) for f in files}
            self.assertIn('feature-summary.tsv', names)
            summary = next(os.path.join(r, f) for r, _, fs in os.walk(d)
                           for f in fs if f == 'feature-summary.tsv')
            df = pd.read_csv(summary, sep='\t', index_col=0)
            self.assertEqual(sorted(df.columns), [
                'direction', 'importance', 'log2-metagenome-over-amplicon',
                'mean-amplicon', 'mean-metagenome', 'prevalence-amplicon',
                'prevalence-metagenome', 'rank'])
            self.assertEqual(set(df['direction']), {'shared'})
            self.assertEqual(df['prevalence-amplicon'].tolist(), [6, 6, 6])

    def test_unpaired_samples_dropped(self):
        md = self.metadata.view(qiime2.Metadata).to_dataframe()
        md = md[md.index != 'sample-6-metagenome']   # sample-6 loses its pair
        metadata = qiime2.Artifact.import_data('ImmutableMetadata',
                                               qiime2.Metadata(md))
        _, dms, _ = self.plugin.pipelines['paired_report'](
            self.table, self.tree, metadata, n_boots=2, permutations=99,
            cv=3, n_estimators=50, random_state=1)
        dm = dms['braycurtis'].view(__import__('skbio').DistanceMatrix)
        self.assertEqual(len(dm.ids), 10)
        self.assertNotIn('sample-6-amplicon', dm.ids)

    def test_too_few_pairs(self):
        md = self.metadata.view(qiime2.Metadata).to_dataframe()
        keep = ~md.index.str.endswith('-metagenome')
        md = md[keep | (md.index == 'sample-1-metagenome')]
        metadata = qiime2.Artifact.import_data('ImmutableMetadata',
                                               qiime2.Metadata(md))
        with self.assertRaisesRegex(ValueError, 'at least 2 are needed'):
            self.plugin.pipelines['paired_report'](
                self.table, self.tree, metadata, n_boots=2, permutations=99,
                cv=3, n_estimators=50)

    def test_wrong_metadata(self):
        md = ex.paired_metadata_factory()
        with self.assertRaisesRegex(ValueError, 'no .original-sample-id. column'):
            self.plugin.pipelines['paired_report'](self.table, self.tree, md)
