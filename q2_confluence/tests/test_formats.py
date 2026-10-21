# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest

import skbio
from qiime2 import Artifact
from qiime2.plugin import ValidationError
from q2_types.tree import NewickFormat

from q2_confluence._formats import GTDBNewickFormat, _drop_root_label

# GTDB's shape: bootstrap:taxon internal labels, a quoted root label with no
# branch length, a space inside a quoted label, a bare bootstrap.
GTDB = ("((RS_GCF_1:0.1,RS_GCF_2:0.3)'100:f__A; g__B':0.5,GB_GCA_3:0.4)"
        "'100:d__Bacteria';\n")


class DropRootLabelTests(unittest.TestCase):

    def test_quoted_root_label(self):
        # TODO: this is AI-generated code; human review pending
        self.assertEqual(
            _drop_root_label(GTDB),
            "((RS_GCF_1:0.1,RS_GCF_2:0.3)'100:f__A; g__B':0.5,GB_GCA_3:0.4);\n")

    def test_bare_root_label(self):
        # TODO: this is AI-generated code; human review pending
        self.assertEqual(_drop_root_label("(a:0.1,b:0.2)100;"),
                         "(a:0.1,b:0.2);\n")

    def test_root_with_branch_length_kept(self):
        # TODO: this is AI-generated code; human review pending
        tree = "(a:0.1,b:0.2)'100:d__Bacteria':0.0;"
        self.assertEqual(_drop_root_label(tree), tree + '\n')

    def test_no_root_label_unchanged(self):
        # TODO: this is AI-generated code; human review pending
        self.assertEqual(_drop_root_label("(a:0.1,b:0.2);"), "(a:0.1,b:0.2);\n")


class GTDBNewickFormatTests(unittest.TestCase):

    def setUp(self):
        # TODO: this is AI-generated code; human review pending
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, 'bac120.tree')
        with open(self.path, 'w') as fh:
            fh.write(GTDB)

    def tearDown(self):
        # TODO: this is AI-generated code; human review pending
        self.tmp.cleanup()

    def test_validates(self):
        # TODO: this is AI-generated code; human review pending
        GTDBNewickFormat(self.path, mode='r').validate('max')

    def test_rejects_non_newick(self):
        # TODO: this is AI-generated code; human review pending
        with open(self.path, 'w') as fh:
            fh.write('Feature ID\tTaxon\n')
        with self.assertRaisesRegex(ValidationError, 'does not start'):
            GTDBNewickFormat(self.path, mode='r').validate('min')
        with open(self.path, 'w') as fh:
            fh.write('(a:0.1,b:0.2')
        with self.assertRaisesRegex(ValidationError, 'does not end'):
            GTDBNewickFormat(self.path, mode='r').validate('max')

    def test_import_drops_the_root_label_only(self):
        # TODO: this is AI-generated code; human review pending
        tree = Artifact.import_data('Phylogeny[Rooted]', self.path,
                                    GTDBNewickFormat)
        newick = tree.view(NewickFormat)
        with newick.open() as fh:
            text = fh.read()
        self.assertTrue(text.rstrip().endswith('GB_GCA_3:0.4);'), text)
        parsed = skbio.TreeNode.read(str(newick), format='newick',
                                     convert_underscores=False)
        self.assertEqual({t.name for t in parsed.tips()},
                         {'RS_GCF_1', 'RS_GCF_2', 'GB_GCA_3'})
        self.assertIsNone(parsed.name)   # the root label is gone
        self.assertEqual({n.name for n in parsed.non_tips()},
                         {'100:f__A; g__B'})   # the internal label stays


if __name__ == '__main__':
    unittest.main()
