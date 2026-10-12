# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import qiime2

from q2_confluence._metadata import merge_metadata


def _md(rows, **cols):
    df = pd.DataFrame(cols, index=pd.Index(rows, name='sample-id'))
    return qiime2.Metadata(df)


class MergeMetadataTests(unittest.TestCase):

    def setUp(self):
        self.a = _md(['s1', 's2', 's3'], site=['x', 'y', 'z'], ph=[6.5, 7.0, 7.5])
        self.b = _md(['s2', 's3', 's4'], site=['y', 'z', 'w'], week=[1, 2, 3])

    def test_union(self):
        out = merge_metadata({'a': self.a, 'b': self.b}).to_dataframe()
        self.assertEqual(list(out.index), ['s1', 's2', 's3', 's4'])
        self.assertEqual(sorted(out.columns), ['ph', 'site', 'week'])
        self.assertEqual(out.loc['s4', 'site'], 'w')
        self.assertTrue(pd.isna(out.loc['s1', 'week']))
        self.assertEqual(out.loc['s2', 'week'], 1)

    def test_intersection(self):
        out = merge_metadata({'a': self.a, 'b': self.b}, 'intersection').to_dataframe()
        self.assertEqual(list(out.index), ['s2', 's3'])
        self.assertEqual(sorted(out.columns), ['ph', 'site', 'week'])

    def test_empty_intersection(self):
        c = _md(['s9'], site=['q'])
        with self.assertRaisesRegex(ValueError, 'nothing to intersect'):
            merge_metadata({'a': self.a, 'c': c}, 'intersection')

    def test_conflict(self):
        b = _md(['s2'], site=['other'])
        with self.assertRaisesRegex(ValueError, "column 'site' has different values"):
            merge_metadata({'a': self.a, 'b': b})

    def test_single_input(self):
        out = merge_metadata({'a': self.a}).to_dataframe()
        pd.testing.assert_frame_equal(out, self.a.to_dataframe())
