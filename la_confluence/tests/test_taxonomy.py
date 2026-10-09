# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from la_confluence._taxonomy import collapsed_name, is_named


class CollapsedNameTests(unittest.TestCase):

    def test_matches_q2_taxa_rule(self):
        taxon = 'd__Bacteria; p__Bacillota; c__Bacilli'
        self.assertEqual(collapsed_name(taxon, 2, 7), 'd__Bacteria;p__Bacillota')
        # padded with '__' to the taxonomy's rank count before cutting
        self.assertEqual(collapsed_name(taxon, 5, 7),
                         'd__Bacteria;p__Bacillota;c__Bacilli;__;__')
        self.assertEqual(collapsed_name(taxon, 3, 3),
                         'd__Bacteria;p__Bacillota;c__Bacilli')


class IsNamedTests(unittest.TestCase):

    def test_named(self):
        self.assertTrue(is_named('d__Bacteria;p__Bacillota'))
        self.assertTrue(is_named('g__Blautia_A'))

    def test_unnamed(self):
        for name in ['d__Bacteria;__', 'd__Bacteria;g__', 'd__Bacteria;',
                     '__', 'g__', '']:
            self.assertFalse(is_named(name), name)
