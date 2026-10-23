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

from qiime2 import Artifact
from q2_types.tree import NewickFormat

from q2_confluence._formats import GTDBNewickFormat

from q2_confluence._tree import _collapse_tree, gtdb_tree_collapse


def _tree(newick):
    # TODO: this is AI-generated code; human review pending
    return skbio.TreeNode.read([newick], convert_underscores=False)


class CollapseTreeTests(unittest.TestCase):

    def setUp(self):
        # GTDB-style: bootstrap:taxon labels on internal nodes, genome tips.
        #                               /-g1
        #                     /100.0:g__A
        #           /90.0:f__F          \-g2
        #          |         |
        # -root----|          \-g3
        #          |
        #          |          /-g4
        #           \95.0:g__C
        #                     \-g5
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
        # Tips named by the collapsed string, shown here by their last rank:
        #                     /-g__A
        #           /90.0:f__F
        # -root----|          \-g__B
        #          |
        #           \-g__C
        out, log = _collapse_tree(self.tree, self.taxonomy, 6)
        tips = {t.name: t for t in out.tips()}
        self.assertEqual(set(tips), {'d__B;p__P;c__C;o__O;f__F;g__A',
                                     'd__B;p__P;c__C;o__O;f__F;g__B',
                                     'd__B;p__P;c__C;o__O;f__G;g__C'})
        # the ancestor keeps its own branch length; nothing is added for
        # the tips sheared off below it
        self.assertAlmostEqual(tips['d__B;p__P;c__C;o__O;f__F;g__A'].length, 0.5)
        # a single-genome genus keeps its branch
        self.assertAlmostEqual(tips['d__B;p__P;c__C;o__O;f__F;g__B'].length, 0.4)
        self.assertAlmostEqual(tips['d__B;p__P;c__C;o__O;f__G;g__C'].length, 2.0)
        # the other internal labels are kept
        self.assertEqual({n.name for n in out.non_tips(include_self=True)},
                         {'90.0:f__F', 'root'})
        # the log: one row per tip, keyed by the new name
        log = log.to_dataframe()
        self.assertEqual(list(log.index), sorted(tips))
        self.assertEqual(log.index.name, 'feature-id')
        self.assertEqual(list(log.columns),
                         ['original-name', 'genomes', 'genome-ids'])
        sheared = log.loc['d__B;p__P;c__C;o__O;f__F;g__A']
        self.assertEqual(sheared['original-name'], '100.0:g__A')
        self.assertEqual(sheared['genomes'], 2)
        self.assertEqual(sheared['genome-ids'], 'g1,g2')
        renamed = log.loc['d__B;p__P;c__C;o__O;f__F;g__B']
        self.assertEqual(renamed['original-name'], 'g3')
        self.assertEqual(renamed['genomes'], 1)
        self.assertEqual(renamed['genome-ids'], 'g3')

    def test_species_is_a_rename(self):
        # TODO: this is AI-generated code; human review pending
        out, log = _collapse_tree(self.tree, self.taxonomy, 7)
        self.assertEqual(out.count(tips=True), 5)
        self.assertIn('d__B;p__P;c__C;o__O;f__F;g__A;s__A one',
                      {t.name for t in out.tips()})

    def test_family_and_root_collapse(self):
        # TODO: this is AI-generated code; human review pending
        #           /-f__F
        # -root----|
        #           \-f__G
        out, log = _collapse_tree(self.tree, self.taxonomy, 5)
        self.assertEqual({t.name for t in out.tips()},
                         {'d__B;p__P;c__C;o__O;f__F', 'd__B;p__P;c__C;o__O;f__G'})
        # a level that leaves one taxon collapses the root onto it
        out, log = _collapse_tree(self.tree, self.taxonomy, 1)
        self.assertTrue(out.is_tip())
        self.assertEqual(out.name, 'd__B')
        log = log.to_dataframe()
        self.assertEqual(log.loc['d__B', 'original-name'], 'root')
        self.assertEqual(log.loc['d__B', 'genomes'], 5)
        self.assertEqual(log.loc['d__B', 'genome-ids'], 'g1,g2,g3,g4,g5')

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
        #                     /-RS_GCF_1
        #           /100.0:g__A
        # ---------|          \-RS_GCF_2
        #          |
        #           \-GB_GCA_3
        # imported as a GTDB tree is, so the tips are quoted and the stock
        # skbio.TreeNode view, which the action takes, keeps their names
        fmt = GTDBNewickFormat()
        with fmt.open() as fh:
            fh.write("((RS_GCF_1:0.1,RS_GCF_2:0.3)'100.0:g__A':0.5,GB_GCA_3:0.4);")
        tree = Artifact.import_data('Phylogeny[Rooted]', fmt).view(skbio.TreeNode)
        lineage = 'd__B;p__P;c__C;o__O;f__F;'
        taxonomy = pd.Series({'RS_GCF_1': lineage + 'g__Blautia_A;s__x',
                              'RS_GCF_2': lineage + 'g__Blautia_A;s__y',
                              'GB_GCA_3': lineage + 'g__Other;s__z'})
        out, log = gtdb_tree_collapse(tree, taxonomy, 6)
        self.assertEqual({t.name for t in out.tips()},
                         {lineage + 'g__Blautia_A', lineage + 'g__Other'})
        # written quoted, so q2-types' default reader keeps the underscores
        fmt2 = NewickFormat()
        with fmt2.open() as fh:
            out.write(fh, format='newick')
        back = skbio.TreeNode.read(str(fmt2), format='newick')
        self.assertEqual({t.name for t in back.tips()},
                         {lineage + 'g__Blautia_A', lineage + 'g__Other'})
        # the kept internal label survives the round trip as well
        out2, _ = gtdb_tree_collapse(tree, taxonomy, 7)
        fmt3 = NewickFormat()
        with fmt3.open() as fh:
            out2.write(fh, format='newick')
        back = skbio.TreeNode.read(str(fmt3), format='newick')
        self.assertEqual({n.name for n in back.non_tips()}, {'100.0:g__A'})
