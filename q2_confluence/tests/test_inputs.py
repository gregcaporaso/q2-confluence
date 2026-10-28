# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import tempfile
import unittest

import biom
import numpy as np
import pandas as pd
import qiime2
import skbio

from q2_confluence._inputs import summarize_inputs


class SummarizeInputsTests(unittest.TestCase):

    def setUp(self):
        # TODO: this is AI-generated code; human review pending
        self.t16s = biom.Table(np.ones((2, 3)), ['f1', 'f2'], ['s1', 's2', 's3'])
        self.tmg = biom.Table(np.ones((3, 2)), ['f2', 'f3', 'f4'], ['s2', 's3'])
        self.tree = skbio.TreeNode.read(io.StringIO('((f1,f2),f3);'))
        self.md = qiime2.Metadata(pd.DataFrame(
            {'site': ['x', 'y', 'z', 'w']},
            index=pd.Index(['s1', 's2', 's3', 's4'], name='sample-id')))

    def test_page_and_membership(self):
        # TODO: this is AI-generated code; human review pending
        with tempfile.TemporaryDirectory() as out:
            summarize_inputs(out, {'16s': self.t16s, 'mg': self.tmg}, self.tree,
                             {'16s': self.md})
            html = open(os.path.join(out, 'index.html')).read()
            self.assertIn('4 sample IDs in all; 2 present in every table', html)
            # f2 is in both tables and the tree; f4 is in no tree
            self.assertIn('4 feature IDs in all; 1 present in every table '
                          'and the tree', html)
            self.assertIn('<td>0.6667</td>', html)   # mg: f2, f3 of f2, f3, f4
            m = pd.read_csv(os.path.join(out, 'sample-membership.tsv'),
                            sep='\t', index_col=0)
            self.assertEqual(list(m.index), ['s1', 's2', 's3', 's4'])
            self.assertEqual(list(m.columns),
                             ['table-16s', 'table-mg', 'metadata-16s'])
            self.assertEqual(m.loc['s4'].tolist(), [False, False, True])
            self.assertEqual(m.loc['s2'].tolist(), [True, True, True])

    def test_tables_only(self):
        # TODO: this is AI-generated code; human review pending
        with tempfile.TemporaryDirectory() as out:
            summarize_inputs(out, {'16s': self.t16s, 'mg': self.tmg}, self.tree)
            html = open(os.path.join(out, 'index.html')).read()
            self.assertIn('3 sample IDs in all; 2 present in every table', html)

    def test_metadata_keyed_like_the_tables(self):
        # TODO: this is AI-generated code; human review pending
        with tempfile.TemporaryDirectory() as out:
            with self.assertRaisesRegex(ValueError, r"metadata keys \['study'\]"):
                summarize_inputs(out, {'16s': self.t16s}, self.tree,
                                 {'study': self.md})
