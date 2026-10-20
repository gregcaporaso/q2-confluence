# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import skbio

from q2_types.tree import NewickFormat

from q2_confluence._tree import _collapse_tree, gtdb_tree_collapse


def _tree(newick):
    # TODO: this is AI-generated code; human review pending
    return skbio.TreeNode.read([newick], convert_underscores=False)


class CollapseTreeTests(unittest.TestCase):

    def setUp(self):
        # GTDB-style: bootstrap:taxon labels on internal nodes, genome tips.
        # TODO: this is AI-generated code; human review pending
        self.tree = _tree(
            "(((g1:0.1,g2:0.3)'100.0:g__A':0.5,g3:0.4)'90.0:f__F':1.0,"
            "(g4:0.2,g5:0.2)'95.0:g__C':2.0)root;")
        self.taxonomy = pd.Series({
            'g1': 'd__B;p__P;c__C;o__O;f__F;g__A;s__A one',
            'g2': 'd__B;p__P;c__C;o__O;f__F;g__A;s__A two',
            'g3': 'd__B;p__P;c__C;o__O;f__F;g__B;s__B one',
            'g4': 'd__B;p__P;c__C;o__O;f__G;g__C;s__C one',
            'g5': 'd__B;p__P;c__C;o__O;f__G;g__C;s__C two',
        })

    def test_genus(self):
        # TODO: this is AI-generated code; human review pending
        out = _collapse_tree(self.tree, self.taxonomy, 6)
        tips = {t.name: t for t in out.tips()}
        self.assertEqual(set(tips), {'d__B;p__P;c__C;o__O;f__F;g__A',
                                     'd__B;p__P;c__C;o__O;f__F;g__B',
                                     'd__B;p__P;c__C;o__O;f__G;g__C'})
        # stem 0.5 plus the median tip depth (0.1 + 0.3) / 2
        self.assertAlmostEqual(tips['d__B;p__P;c__C;o__O;f__F;g__A'].length, 0.7)
        # a single-genome genus keeps its branch
        self.assertAlmostEqual(tips['d__B;p__P;c__C;o__O;f__F;g__B'].length, 0.4)
        self.assertAlmostEqual(tips['d__B;p__P;c__C;o__O;f__G;g__C'].length, 2.2)
        self.assertTrue(all(n.name is None for n in out.non_tips(include_self=True)))

    def test_species_is_a_rename(self):
        # TODO: this is AI-generated code; human review pending
        out = _collapse_tree(self.tree, self.taxonomy, 7)
        self.assertEqual(out.count(tips=True), 5)
        self.assertIn('d__B;p__P;c__C;o__O;f__F;g__A;s__A one',
                      {t.name for t in out.tips()})

    def test_family_and_root_collapse(self):
        # TODO: this is AI-generated code; human review pending
        out = _collapse_tree(self.tree, self.taxonomy, 5)
        self.assertEqual({t.name for t in out.tips()},
                         {'d__B;p__P;c__C;o__O;f__F', 'd__B;p__P;c__C;o__O;f__G'})
        # a level that leaves one taxon collapses the root onto it
        out = _collapse_tree(self.tree, self.taxonomy, 1)
        self.assertTrue(out.is_tip())
        self.assertEqual(out.name, 'd__B')

    def test_input_tree_unchanged(self):
        # TODO: this is AI-generated code; human review pending
        before = str(self.tree)
        _collapse_tree(self.tree, self.taxonomy, 6)
        self.assertEqual(str(self.tree), before)

    def test_not_monophyletic(self):
        # g2 and g3 as one genus: their common ancestor also holds g1
        # TODO: this is AI-generated code; human review pending
        taxonomy = self.taxonomy.copy()
        taxonomy['g2'] = 'd__B;p__P;c__C;o__O;f__F;g__B;s__B two'
        with self.assertRaisesRegex(ValueError, 'not monophyletic'):
            _collapse_tree(self.tree, taxonomy, 6)

    def test_missing_tip(self):
        # TODO: this is AI-generated code; human review pending
        with self.assertRaisesRegex(ValueError, 'missing from the taxonomy'):
            _collapse_tree(self.tree, self.taxonomy.drop('g5'), 6)

    def test_unnamed_at_level(self):
        # TODO: this is AI-generated code; human review pending
        taxonomy = self.taxonomy.copy()
        taxonomy['g3'] = 'd__B;p__P;c__C;o__O;f__F;g__;s__'
        with self.assertRaisesRegex(ValueError, 'no name at level 6'):
            _collapse_tree(self.tree, taxonomy, 6)

    def test_bad_level(self):
        # TODO: this is AI-generated code; human review pending
        for level in [0, 8]:
            with self.assertRaisesRegex(ValueError, 'between 1 and 7'):
                _collapse_tree(self.tree, self.taxonomy, level)

    def test_not_seven_ranks(self):
        # TODO: this is AI-generated code; human review pending
        taxonomy = self.taxonomy.copy()
        taxonomy['g3'] = 'd__B;p__P;c__C;o__O;f__F;g__B'
        with self.assertRaisesRegex(ValueError, '1 taxonomy strings do not '
                                    'have the 7 GTDB ranks'):
            _collapse_tree(self.tree, taxonomy, 6)


class GtdbTreeCollapseTests(unittest.TestCase):

    def test_underscores_kept_and_names_round_trip(self):
        # TODO: this is AI-generated code; human review pending
        fmt = NewickFormat()
        with fmt.open() as fh:
            fh.write("((RS_GCF_1:0.1,RS_GCF_2:0.3)'100.0:g__A':0.5,GB_GCA_3:0.4);")
        lineage = 'd__B;p__P;c__C;o__O;f__F;'
        taxonomy = pd.Series({'RS_GCF_1': lineage + 'g__Blautia_A;s__x',
                              'RS_GCF_2': lineage + 'g__Blautia_A;s__y',
                              'GB_GCA_3': lineage + 'g__Other;s__z'})
        out = gtdb_tree_collapse(fmt, taxonomy, 6)
        self.assertEqual({t.name for t in out.tips()},
                         {lineage + 'g__Blautia_A', lineage + 'g__Other'})
        # written quoted, so q2-types' default reader keeps the underscores
        fmt2 = NewickFormat()
        with fmt2.open() as fh:
            out.write(fh, format='newick')
        back = skbio.TreeNode.read(str(fmt2), format='newick')
        self.assertEqual({t.name for t in back.tips()},
                         {lineage + 'g__Blautia_A', lineage + 'g__Other'})
