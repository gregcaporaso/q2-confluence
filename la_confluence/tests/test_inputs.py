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

import biom
import numpy as np
import pandas as pd
import qiime2

from la_confluence._inputs import summarize_inputs


class SummarizeInputsTests(unittest.TestCase):

    def setUp(self):
        self.t16s = biom.Table(np.ones((2, 3)), ['f1', 'f2'], ['s1', 's2', 's3'])
        self.tmg = biom.Table(np.ones((3, 2)), ['g1', 'g2', 'g3'], ['s2', 's3'])
        self.md = qiime2.Metadata(pd.DataFrame(
            {'site': ['x', 'y', 'z', 'w']},
            index=pd.Index(['s1', 's2', 's3', 's4'], name='sample-id')))

    def test_page_and_membership(self):
        with tempfile.TemporaryDirectory() as out:
            summarize_inputs(out, {'16s': self.t16s, 'mg': self.tmg},
                             {'study': self.md})
            html = open(os.path.join(out, 'index.html')).read()
            self.assertIn('4 sample IDs in all; 2 present in every input', html)
            m = pd.read_csv(os.path.join(out, 'sample-membership.tsv'),
                            sep='\t', index_col=0)
            self.assertEqual(list(m.index), ['s1', 's2', 's3', 's4'])
            self.assertEqual(list(m.columns),
                             ['table 16s', 'table mg', 'metadata study'])
            self.assertEqual(m.loc['s4'].tolist(), [False, False, True])
            self.assertEqual(m.loc['s2'].tolist(), [True, True, True])

    def test_tables_only(self):
        with tempfile.TemporaryDirectory() as out:
            summarize_inputs(out, {'16s': self.t16s, 'mg': self.tmg})
            html = open(os.path.join(out, 'index.html')).read()
            self.assertIn('3 sample IDs in all; 2 present in every input', html)
