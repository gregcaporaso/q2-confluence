# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile
import unittest

import numpy as np

from q2_confluence._sylph import (_accession, _profiles_to_artifacts,
                                  _query_to_metadata, _read_sylphmpa)

LINEAGE = ('d__Bacteria|p__Bacillota_A|c__Clostridia|o__Lachnospirales|'
           'f__Lachnospiraceae|g__Blautia')
HEADER = ('clade_name\trelative_abundance\tsequence_abundance\t'
          'ANI (if strain-level)\tCoverage (if strain-level)\n')


def _sylphmpa(sample, strains):
    # TODO: this is AI-generated code; human review pending
    """strains: [(accession, species, rel, seq, ani)] -> file text."""
    lines = [f'#SampleID\t{sample}\tTaxonomies_used:[\'GTDB_r232\']\n', HEADER,
             f'd__Bacteria\t100.0\t{sum(s[3] for s in strains)}\tNA\tNA\n']
    for acc, species, rel, seq, ani in strains:
        lines.append(f'{LINEAGE}|s__{species}\t{rel}\t{seq}\tNA\tNA\n')
        lines.append(f'{LINEAGE}|s__{species}|t__{acc}_genomic.fna.gz\t{rel}\t{seq}'
                     f'\t{ani}\t12.5\n')
    return ''.join(lines)


class AccessionTests(unittest.TestCase):

    def test_forms(self):
        # TODO: this is AI-generated code; human review pending
        self.assertEqual(_accession('GCF_000227485.1_ASM22748v1_genomic.fna.gz'),
                         'GCF_000227485.1')
        self.assertEqual(_accession('/db/GCA_01.2_genomic.fna'), 'GCA_01.2')
        self.assertEqual(_accession('mygenome.fasta'), 'mygenome')


class ProfilesTests(unittest.TestCase):

    def setUp(self):
        # TODO: this is AI-generated code; human review pending
        self.tmp = tempfile.TemporaryDirectory()
        self.paths = []
        for sample, strains in [
            ('s1', [('GCF_1.1', 'Blautia wexlerae', 60.0, 70.0, 99.1),
                    ('GCF_2.1', 'Blautia_A obeum', 40.0, 10.0, 97.3)]),
            ('s2', [('GCF_2.1', 'Blautia_A obeum', 100.0, 50.0, 98.0)]),
        ]:
            p = os.path.join(self.tmp.name, f'{sample}.sylphmpa')
            open(p, 'w').write(_sylphmpa(sample, strains))
            self.paths.append(p)

    def tearDown(self):
        # TODO: this is AI-generated code; human review pending
        self.tmp.cleanup()

    def test_read_one(self):
        # TODO: this is AI-generated code; human review pending
        sample, df = _read_sylphmpa(self.paths[0])
        self.assertEqual(sample, 's1')
        self.assertEqual(list(df.index), ['GCF_1.1', 'GCF_2.1'])
        self.assertEqual(df.loc['GCF_1.1', 'Taxon'],
                         LINEAGE.replace('|', ';') + ';s__Blautia wexlerae')
        self.assertEqual(df.loc['GCF_2.1', 'ani'], 97.3)

    def test_artifacts(self):
        # TODO: this is AI-generated code; human review pending
        rel, seq, tax, stats = _profiles_to_artifacts(self.paths)
        self.assertEqual(list(rel.ids()), ['s1', 's2'])
        self.assertEqual(list(rel.ids(axis='observation')), ['GCF_1.1', 'GCF_2.1'])
        self.assertAlmostEqual(rel.get_value_by_ids('GCF_1.1', 's1'), 0.6)
        self.assertAlmostEqual(rel.get_value_by_ids('GCF_1.1', 's2'), 0.0)
        # sequence abundance renormalized over detected genomes (unknown 20%)
        self.assertAlmostEqual(seq.get_value_by_ids('GCF_1.1', 's1'), 0.875)
        self.assertAlmostEqual(seq.get_value_by_ids('GCF_2.1', 's2'), 1.0)
        self.assertEqual(tax.index.name, 'Feature ID')
        self.assertEqual(tax['GCF_2.1'],
                         LINEAGE.replace('|', ';') + ';s__Blautia_A obeum')
        self.assertEqual(stats.loc['s1', 'genomes'], 2)
        self.assertEqual(stats.loc['s1', 'genera'], 1)
        self.assertAlmostEqual(stats.loc['s1', 'unknown-share'], 0.2)
        self.assertAlmostEqual(stats.loc['s2', 'unknown-share'], 0.5)
        self.assertAlmostEqual(stats.loc['s1', 'median-ani'], 98.2)

    def test_conflicting_taxonomy(self):
        # TODO: this is AI-generated code; human review pending
        p = os.path.join(self.tmp.name, 's3.sylphmpa')
        open(p, 'w').write(_sylphmpa('s3', [('GCF_1.1', 'Other name', 100, 50, 99)]))
        with self.assertRaisesRegex(ValueError, 'two taxonomy strings'):
            _profiles_to_artifacts(self.paths + [p])

    def test_not_sylphmpa(self):
        # TODO: this is AI-generated code; human review pending
        p = os.path.join(self.tmp.name, 'x.tsv')
        open(p, 'w').write('a\tb\n')
        with self.assertRaisesRegex(ValueError, 'not a .sylphmpa'):
            _read_sylphmpa(p)


_QUERY_HEADER = ('Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\t'
                 'ANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\t'
                 'Median_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\t'
                 'Contig_name\n')


class QueryParserTests(unittest.TestCase):

    def setUp(self):
        # TODO: this is AI-generated code; human review pending
        self.dir = tempfile.mkdtemp()
        self.path = os.path.join(self.dir, 'query.tsv')
        with open(self.path, 'w') as fh:
            fh.write(_QUERY_HEADER)
            fh.write('s2\t/db/GCF_000002.1_genomic.fna.gz\t97.50\t0.080\t96.1-101.2\t'
                     '0.08\t0.05-0.11\t1\t1.000\t60/900\t95.10\tGCF_000002.1\n')
            fh.write('s1\t/db/GCF_000001.1_genomic.fna.gz\t100.00\t15.000\tNA-NA\t'
                     'HIGH\tNA-NA\t15\t15.485\t270/270\t100.00\tGCF_000001.1\n')
        self.taxonomy = {'GCF_000001.1': 'd__Bacteria;g__A;s__A one',
                         'GCF_000002.1': 'd__Bacteria;g__B;s__B two'}

    def tearDown(self):
        # TODO: this is AI-generated code; human review pending
        shutil.rmtree(self.dir)

    def test_rows(self):
        # TODO: this is AI-generated code; human review pending
        df = _query_to_metadata(self.path, self.taxonomy).to_dataframe()
        self.assertEqual(list(df.index), ['s1-GCF_000001.1', 's2-GCF_000002.1'])
        self.assertEqual(list(df.columns), [
            'sample', 'accession', 'taxonomy', 'adjusted-ani',
            'effective-coverage', 'ani-5-percentile', 'ani-95-percentile',
            'containment-index', 'naive-ani'])
        self.assertEqual(df.loc['s2-GCF_000002.1', 'taxonomy'],
                         'd__Bacteria;g__B;s__B two')
        self.assertEqual(df.loc['s2-GCF_000002.1', 'adjusted-ani'], 97.5)
        self.assertEqual(df.loc['s2-GCF_000002.1', 'effective-coverage'], 0.08)
        self.assertEqual(df.loc['s2-GCF_000002.1', 'ani-5-percentile'], 96.1)
        self.assertEqual(df.loc['s2-GCF_000002.1', 'ani-95-percentile'], 101.2)
        self.assertTrue(np.isnan(df.loc['s1-GCF_000001.1', 'ani-5-percentile']))
        self.assertEqual(df.loc['s1-GCF_000001.1', 'containment-index'], '270/270')
        self.assertEqual(df.loc['s1-GCF_000001.1', 'naive-ani'], 100.0)

    def test_targets(self):
        # TODO: this is AI-generated code; human review pending
        by_taxon = _query_to_metadata(self.path, self.taxonomy, ['g__B;'])
        self.assertEqual(list(by_taxon.to_dataframe().index), ['s2-GCF_000002.1'])
        by_accession = _query_to_metadata(self.path, self.taxonomy, ['GCF_000001.1'])
        self.assertEqual(list(by_accession.to_dataframe().index), ['s1-GCF_000001.1'])
        with self.assertRaisesRegex(ValueError, 'no reported genome matches'):
            _query_to_metadata(self.path, self.taxonomy, ['g__C;'])

    def test_empty(self):
        # TODO: this is AI-generated code; human review pending
        with open(self.path, 'w') as fh:
            fh.write(_QUERY_HEADER)
        with self.assertRaisesRegex(ValueError, 'no genome with enough'):
            _query_to_metadata(self.path, self.taxonomy)
