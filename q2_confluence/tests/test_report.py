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
        # TODO: this is AI-generated code; human review pending
        super().setUp()
        # the genome space: the reference tree as it is, no level
        self.tree = ex._reference_tree_factory()
        # the reference taxonomy is used as it is: integrate no longer
        # produces one, and the report takes any taxonomy covering the table
        self.taxonomy = ex._reference_taxonomy_factory()
        self.table, self.metadata = self.plugin.pipelines['integrate'](
            {'amplicon': ex._paired_amplicon_table_factory(),
             'metagenome': ex._paired_metagenome_table_factory()},
            self.tree, metadata={'amplicon': ex._paired_metadata_factory(),
                                 'metagenome': ex._paired_metadata_factory()})

    def test_report(self):
        # TODO: this is AI-generated code; human review pending
        report, dms, importance, summary = self.plugin.pipelines['paired_report'](
            self.table, self.tree, self.metadata, taxonomy=self.taxonomy,
            permutations=99, cv=3, n_estimators=50, random_state=1)
        self.assertEqual(sorted(dms.keys()),
                         ['braycurtis', 'jaccard', 'unweighted_unifrac',
                          'weighted_unifrac'])
        imp = importance.view(pd.DataFrame)
        self.assertEqual(len(imp), 3)
        df = summary.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(df.index.name, 'feature-id')
        self.assertEqual(sorted(df.columns), [
            'direction', 'importance', 'log2-metagenome-over-amplicon',
            'mean-amplicon', 'mean-metagenome', 'prevalence-amplicon',
            'prevalence-metagenome', 'rank'])
        self.assertEqual(set(df['direction']), {'shared'})
        self.assertEqual(df['prevalence-amplicon'].tolist(), [6, 6, 6])
        self.assertEqual(df['rank'].tolist(), [1, 2, 3])
        with tempfile.TemporaryDirectory() as d:
            report.export_data(d)
            tabulated = [os.path.join(root, f) for root, _, files in os.walk(d)
                         for f in files if f == 'metadata.tsv']
            self.assertEqual(len(tabulated), 1)   # the tabulated summary
            with open(tabulated[0]) as fh:
                header = fh.readline().rstrip('\n').split('\t')
            self.assertEqual(header[-1], 'Taxon')   # the taxonomy joined

    def test_report_without_taxonomy(self):
        # TODO: this is AI-generated code; human review pending
        report, _, _, _ = self.plugin.pipelines['paired_report'](
            self.table, self.tree, self.metadata, permutations=99,
            cv=3, n_estimators=50, random_state=1)
        with tempfile.TemporaryDirectory() as d:
            report.export_data(d)
            tabulated = [os.path.join(root, f) for root, _, files in os.walk(d)
                         for f in files if f == 'metadata.tsv']
            with open(tabulated[0]) as fh:
                header = fh.readline().rstrip('\n').split('\t')
            self.assertNotIn('Taxon', header)

    def test_taxonomy_missing_features(self):
        # TODO: this is AI-generated code; human review pending
        partial = ex._taxonomy_factory({'RS_GCF_1': ex._TAXONOMY['RS_GCF_1']})
        with self.assertRaisesRegex(ValueError, 'not in the taxonomy'):
            self.plugin.pipelines['paired_report'](
                self.table, self.tree, self.metadata, taxonomy=partial,
                permutations=99, cv=3, n_estimators=50)

    def test_summary_directions(self):
        # TODO: this is AI-generated code; human review pending
        summary, = self.plugin.methods['summarize_paired_features'](
            ex._paired_relative_table_factory(), ex._paired_importance_factory(),
            ex._paired_integrated_metadata_factory().view(qiime2.Metadata))
        df = summary.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(df['direction'].tolist(),
                         ['amplicon-only', 'metagenome-only', 'shared'])
        self.assertEqual(df['rank'].tolist(), [1, 2, 3])
        self.assertEqual(df['prevalence-metagenome'].tolist(), [0, 3, 3])
        self.assertAlmostEqual(df['mean-amplicon'].iloc[0], 0.5)

    def test_unpaired_samples_dropped(self):
        # TODO: this is AI-generated code; human review pending
        md = self.metadata.view(qiime2.Metadata).to_dataframe()
        md = md[md.index != 'sample-6-metagenome']   # sample-6 loses its pair
        metadata = qiime2.Artifact.import_data('ImmutableMetadata',
                                               qiime2.Metadata(md))
        _, dms, _, _ = self.plugin.pipelines['paired_report'](
            self.table, self.tree, metadata, permutations=99,
            cv=3, n_estimators=50, random_state=1)
        dm = dms['braycurtis'].view(__import__('skbio').DistanceMatrix)
        self.assertEqual(len(dm.ids), 10)
        self.assertNotIn('sample-6-amplicon', dm.ids)

    def test_too_few_pairs(self):
        # TODO: this is AI-generated code; human review pending
        md = self.metadata.view(qiime2.Metadata).to_dataframe()
        keep = ~md.index.str.endswith('-metagenome')
        md = md[keep | (md.index == 'sample-1-metagenome')]
        metadata = qiime2.Artifact.import_data('ImmutableMetadata',
                                               qiime2.Metadata(md))
        with self.assertRaisesRegex(ValueError, 'at least 2 are needed'):
            self.plugin.pipelines['paired_report'](
                self.table, self.tree, metadata, permutations=99,
                cv=3, n_estimators=50)

    def test_wrong_metadata(self):
        # TODO: this is AI-generated code; human review pending
        md = ex._paired_metadata_factory()
        with self.assertRaisesRegex(ValueError, 'no .original-sample-id. column'):
            self.plugin.pipelines['paired_report'](self.table, self.tree, md)
