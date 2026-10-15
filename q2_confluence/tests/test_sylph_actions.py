# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import os
import tempfile
from unittest import mock

import biom
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_confluence import _examples as ex
from q2_confluence._sylph import SylphDatabaseDirFmt, SylphSketchesDirFmt


class SketchReadsTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def test_sketch(self):
        # TODO: this is AI-generated code; human review pending
        sketches, = self.plugin.methods['sketch_reads'](
            ex._paired_reads_factory(), threads=1)
        self.assertEqual(str(sketches.type), 'SylphSketches')
        paths = sketches.view(SylphSketchesDirFmt).sample_paths()
        self.assertEqual(list(paths), ['sample-1', 'sample-2', 'sample-3'])
        for p in paths.values():
            self.assertTrue(p.endswith('.paired.sylsp'))
            self.assertGreater(os.path.getsize(p), 0)


class ProfileReadsTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def test_profile(self):
        # TODO: this is AI-generated code; human review pending
        sketches, = self.plugin.methods['sketch_reads'](
            ex._paired_reads_factory(), threads=1)
        database = ex._sylph_database_factory()
        rel, seq, tax, stats = self.plugin.methods['profile_reads'](
            sketches, database, threads=1)
        rel = rel.view(biom.Table)
        self.assertEqual(list(rel.ids()), ['sample-1', 'sample-2', 'sample-3'])
        self.assertEqual(list(rel.ids(axis='observation')),
                         ['GCF_000001.1', 'GCF_000002.1'])
        self.assertAlmostEqual(rel.get_value_by_ids('GCF_000001.1', 'sample-1'), 1.0)
        self.assertEqual(rel.get_value_by_ids('GCF_000002.1', 'sample-1'), 0.0)
        self.assertAlmostEqual(rel.get_value_by_ids('GCF_000001.1', 'sample-2'),
                               0.5, delta=0.1)
        # the third sample matches nothing: all zero, unknown share 1
        self.assertEqual(rel.sum(axis='sample')[2], 0.0)
        stats = stats.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(stats.loc['sample-3', 'genomes'], 0)
        self.assertEqual(stats.loc['sample-3', 'unknown-share'], 1.0)
        self.assertEqual(stats.loc['sample-2', 'genomes'], 2)
        self.assertEqual(tax.view(pd.Series)['GCF_000002.1'],
                         ex._TAXONOMY['RS_GCF_1'])
        seq = seq.view(biom.Table)
        self.assertAlmostEqual(seq.sum(axis='sample')[1], 1.0)


class QueryReadsTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def setUp(self):
        # TODO: this is AI-generated code; human review pending
        super().setUp()
        self.sketches, = self.plugin.methods['sketch_reads'](
            ex._paired_reads_factory(), threads=1)
        self.database = ex._sylph_database_factory()
        self.query = self.plugin.methods['query_reads']

    def test_query(self):
        # TODO: this is AI-generated code; human review pending
        query, = self.query(self.sketches, self.database, threads=1)
        df = query.view(qiime2.Metadata).to_dataframe()
        # sample-3 holds reads from a genome the database lacks: no row
        self.assertEqual(list(df.index), ['sample-1-GCF_000001.1',
                                          'sample-2-GCF_000001.1',
                                          'sample-2-GCF_000002.1'])
        self.assertEqual(df.loc['sample-2-GCF_000002.1', 'sample'], 'sample-2')
        self.assertEqual(df.loc['sample-2-GCF_000002.1', 'accession'],
                         'GCF_000002.1')
        self.assertEqual(df.loc['sample-2-GCF_000002.1', 'taxonomy'],
                         ex._TAXONOMY['RS_GCF_1'])
        self.assertAlmostEqual(df.loc['sample-1-GCF_000001.1', 'adjusted-ani'],
                               100.0, delta=0.1)
        self.assertGreater(df.loc['sample-1-GCF_000001.1', 'effective-coverage'],
                           1.0)
        self.assertRegex(df.loc['sample-1-GCF_000001.1', 'containment-index'],
                         r'^\d+/\d+$')

    def test_targets(self):
        # TODO: this is AI-generated code; human review pending
        query, = self.query(self.sketches, self.database, threads=1,
                            targets=['g__Blautia;'])
        df = query.view(qiime2.Metadata).to_dataframe()
        self.assertEqual(list(df.index), ['sample-2-GCF_000002.1'])
        with self.assertRaisesRegex(ValueError, 'no reported genome matches'):
            self.query(self.sketches, self.database, threads=1,
                       targets=['g__Escherichia;'])


class FetchSylphDatabaseTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def test_fetch_from_local_urls(self):
        # TODO: this is AI-generated code; human review pending
        directory = ex._sylph_data()
        db = os.path.join(directory, 'db', 'database.syldb')
        tax = os.path.join(directory, 'db', 'taxonomy.tsv.gz')
        with mock.patch('q2_confluence._sylph._database_url',
                        return_value='file://' + db), \
                mock.patch('q2_confluence._sylph._taxonomy_url',
                           return_value='file://' + tax):
            database, = self.plugin.methods['fetch_sylph_database']('232')
        fmt = database.view(SylphDatabaseDirFmt)
        self.assertTrue(fmt.database_path().endswith('database.syl2db'))
        self.assertEqual(os.path.getsize(fmt.database_path()),
                         os.path.getsize(db))

    def test_real_urls(self):
        # TODO: this is AI-generated code; human review pending
        from q2_confluence._sylph import _database_url, _taxonomy_url
        self.assertEqual(_database_url('232'), 'http://faust.compbio.cs.cmu.edu/'
                         'sylph-stuff/gtdb-r232-c200-dbv2.syl2db')
        self.assertTrue(_taxonomy_url('226').endswith('gtdb_r226_metadata.tsv.gz'))


class SylphDatabaseFormatTests(TestPluginBase):
    package = 'q2_confluence.tests'

    def test_validation(self):
        # TODO: this is AI-generated code; human review pending
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, 'database.syldb'), 'wb').close()
            with gzip.open(os.path.join(d, 'taxonomy.tsv.gz'), 'wt') as fh:
                fh.write('GCF_1\td__Bacteria\n')
            with self.assertRaisesRegex(Exception, 'empty sylph database'):
                SylphDatabaseDirFmt(d, mode='r').validate()

    def test_sketches_validation(self):
        # TODO: this is AI-generated code; human review pending
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, 'sample-1.paired.sylsp'), 'wb').close()
            with self.assertRaisesRegex(Exception, 'empty sylph sketch'):
                SylphSketchesDirFmt(d, mode='r').validate()
