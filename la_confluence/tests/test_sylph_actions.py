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

from la_confluence import _examples as ex
from la_confluence._sylph import SylphDatabaseDirFmt


class ProfileReadsTests(TestPluginBase):
    package = 'la_confluence.tests'

    def test_profile(self):
        reads = ex.paired_reads_factory()
        database = ex.sylph_database_factory()
        rel, seq, tax, stats = self.plugin.methods['profile_reads'](
            reads, database, threads=1)
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


class FetchSylphDatabaseTests(TestPluginBase):
    package = 'la_confluence.tests'

    def test_fetch_from_local_urls(self):
        directory = ex._sylph_data()
        db = os.path.join(directory, 'db', 'database.syldb')
        tax = os.path.join(directory, 'db', 'taxonomy.tsv.gz')
        with mock.patch('la_confluence._sylph._database_url',
                        return_value='file://' + db), \
                mock.patch('la_confluence._sylph._taxonomy_url',
                           return_value='file://' + tax):
            database, = self.plugin.methods['fetch_sylph_database']('232')
        fmt = database.view(SylphDatabaseDirFmt)
        self.assertTrue(fmt.database_path().endswith('database.syl2db'))
        self.assertEqual(os.path.getsize(fmt.database_path()),
                         os.path.getsize(db))

    def test_real_urls(self):
        from la_confluence._sylph import _database_url, _taxonomy_url
        self.assertEqual(_database_url('232'), 'http://faust.compbio.cs.cmu.edu/'
                         'sylph-stuff/gtdb-r232-c200-dbv2.syl2db')
        self.assertTrue(_taxonomy_url('226').endswith('gtdb_r226_metadata.tsv.gz'))


class SylphDatabaseFormatTests(TestPluginBase):
    package = 'la_confluence.tests'

    def test_validation(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, 'database.syldb'), 'wb').close()
            with gzip.open(os.path.join(d, 'taxonomy.tsv.gz'), 'wt') as fh:
                fh.write('GCF_1\td__Bacteria\n')
            with self.assertRaisesRegex(Exception, 'empty sylph database'):
                SylphDatabaseDirFmt(d, mode='r').validate()
