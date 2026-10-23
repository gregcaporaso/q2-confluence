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

from q2_confluence._formats import (GTDBNewickFormat, _drop_root_label,
                                    _quote_underscore_labels)

# GTDB's shape: bootstrap:taxon internal labels, a quoted root label with no
# branch length, a space inside a quoted label, a bare bootstrap.
GTDB = ("((RS_GCF_1:0.1,RS_GCF_2:0.3)'100:f__A; g__B':0.5,GB_GCA_3:0.4)"
        "'100:d__Bacteria';\n")
# The archaeal tree (ar53_r226.tree, ar53_r232.tree) has the same shape; its
# root label is '100:d__Archaea' (checked on both files on 2026-10-22).
ARCHAEA = ("((RS_GCF_4:0.1,GB_GCA_5:0.3)'99.0:p__Halobacteriota':0.5,"
           "GB_GCA_6:0.4)'100:d__Archaea';\n")


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

    def test_archaeal_root_label(self):
        # TODO: this is AI-generated code; human review pending
        self.assertEqual(
            _drop_root_label(ARCHAEA),
            "((RS_GCF_4:0.1,GB_GCA_5:0.3)'99.0:p__Halobacteriota':0.5,"
            "GB_GCA_6:0.4);\n")


class QuoteUnderscoreLabelsTests(unittest.TestCase):

    def test_tips_quoted_once(self):
        # TODO: this is AI-generated code; human review pending
        bare = "((RS_GCF_1:0.1,RS_GCF_2:0.3)100:0.5,GB_GCA_3:0.4);"
        quoted = "(('RS_GCF_1':0.1,'RS_GCF_2':0.3)100:0.5,'GB_GCA_3':0.4);"
        self.assertEqual(_quote_underscore_labels(bare), quoted)
        self.assertEqual(_quote_underscore_labels(quoted), quoted)
        self.assertEqual(_quote_underscore_labels(GTDB),
                         "(('RS_GCF_1':0.1,'RS_GCF_2':0.3)'100:f__A; g__B':0.5,"
                         "'GB_GCA_3':0.4)'100:d__Bacteria';\n")

    def test_labels_without_underscore_untouched(self):
        # TODO: this is AI-generated code; human review pending
        plain = "((a:0.1,b:0.3)100:0.5,c:0.4);"
        self.assertEqual(_quote_underscore_labels(plain), plain)

    def test_internal_label_with_underscore_quoted(self):
        # TODO: this is AI-generated code; human review pending
        self.assertEqual(_quote_underscore_labels("(a:0.1,b:0.3)g__A:0.5;"),
                         "(a:0.1,b:0.3)'g__A':0.5;")


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
        for level in ['min', 'max']:
            GTDBNewickFormat(self.path, mode='r').validate(level)

    def test_rejects_non_newick(self):
        # TODO: this is AI-generated code; human review pending
        # every level parses the whole file, so min rejects what max rejects
        for content, message in [('Feature ID\tTaxon\n', 'does not start'),
                                 ('(a:0.1,b:0.2', 'does not end'),
                                 ('(a:0.1,b:0.2));', 'not a Newick tree')]:
            with open(self.path, 'w') as fh:
                fh.write(content)
            for level in ['min', 'max']:
                with self.assertRaisesRegex(ValidationError, message):
                    GTDBNewickFormat(self.path, mode='r').validate(level)

    def test_import_drops_the_root_label_and_quotes_the_tips(self):
        # TODO: this is AI-generated code; human review pending
        tree = Artifact.import_data('Phylogeny[Rooted]', self.path,
                                    GTDBNewickFormat)
        newick = tree.view(NewickFormat)
        with newick.open() as fh:
            text = fh.read()
        self.assertEqual(
            text, "(('RS_GCF_1':0.1,'RS_GCF_2':0.3)'100:f__A; g__B':0.5,"
                  "'GB_GCA_3':0.4);\n")
        # the stock view (Newick convention) now keeps the underscores
        parsed = tree.view(skbio.TreeNode)
        self.assertEqual({t.name for t in parsed.tips()},
                         {'RS_GCF_1', 'RS_GCF_2', 'GB_GCA_3'})
        self.assertIsNone(parsed.name)   # the root label is gone
        self.assertEqual({n.name for n in parsed.non_tips()},
                         {'100:f__A; g__B'})   # the internal label stays

    def test_import_archaeal_tree(self):
        # TODO: this is AI-generated code; human review pending
        with open(self.path, 'w') as fh:
            fh.write(ARCHAEA)
        tree = Artifact.import_data('Phylogeny[Rooted]', self.path,
                                    GTDBNewickFormat)
        parsed = tree.view(skbio.TreeNode)
        self.assertIsNone(parsed.name)
        self.assertEqual(parsed.count(tips=True), 3)


if __name__ == '__main__':
    unittest.main()
