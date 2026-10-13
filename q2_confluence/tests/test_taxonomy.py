# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

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
