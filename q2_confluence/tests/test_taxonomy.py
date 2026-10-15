# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import biom
import numpy as np
import pandas as pd
from q2_taxa._util import _collapse_table

from q2_confluence._taxonomy import _collapsed_name, _is_named


class CollapsedNameTests(unittest.TestCase):

    def test_matches_q2_taxa_rule(self):
        # TODO: this is AI-generated code; human review pending
        taxon = 'd__Bacteria; p__Bacillota; c__Bacilli'
        self.assertEqual(_collapsed_name(taxon, 2, 7), 'd__Bacteria;p__Bacillota')
        # padded with '__' to the taxonomy's rank count before cutting
        self.assertEqual(_collapsed_name(taxon, 5, 7),
                         'd__Bacteria;p__Bacillota;c__Bacilli;__;__')
        self.assertEqual(_collapsed_name(taxon, 3, 3),
                         'd__Bacteria;p__Bacillota;c__Bacilli')

    def test_equals_q2_taxa_collapse(self):
        # TODO: this is AI-generated code; human review pending
        # the IDs q2-taxa's collapse writes, from its own implementation
        # (q2_taxa._util._collapse_table), for every level of three strings
        taxa = pd.Series({
            'f1': 'd__Bacteria; p__Bacillota; c__Bacilli',
            'f2': 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;'
                  'f__Bacillaceae;g__Bacillus;s__Bacillus subtilis',
            'f3': 'd__Bacteria;p__Bacillota_A;c__;o__;f__;g__;s__'})
        n_ranks = int(taxa.str.split(';').str.len().max())
        table = biom.Table(np.ones((3, 1)), list(taxa.index), ['s1'])
        for level in range(1, n_ranks + 1):
            collapsed = _collapse_table(table, taxa, level, n_ranks)
            self.assertEqual(
                set(collapsed.ids(axis='observation')),
                {_collapsed_name(t, level, n_ranks) for t in taxa}, level)


class IsNamedTests(unittest.TestCase):

    def test_named(self):
        # TODO: this is AI-generated code; human review pending
        self.assertTrue(_is_named('d__Bacteria;p__Bacillota'))
        self.assertTrue(_is_named('g__Blautia_A'))

    def test_unnamed(self):
        # TODO: this is AI-generated code; human review pending
        for name in ['d__Bacteria;__', 'd__Bacteria;g__', 'd__Bacteria;',
                     '__', 'g__', '']:
            self.assertFalse(_is_named(name), name)
