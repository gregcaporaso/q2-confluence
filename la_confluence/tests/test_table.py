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
import numpy.testing as npt

from la_confluence._table import frequency


class FrequencyTests(unittest.TestCase):

    def setUp(self):
        self.table = biom.Table(
            np.array([[0.5, 0.2, 0.0], [0.5, 0.7999996, 0.0], [0.0, 0.0000004, 0.0]]),
            ['f1', 'f2', 'f3'], ['s1', 's2', 'empty'])

    def test_counts(self):
        out = frequency(self.table, 1000)
        npt.assert_array_equal(out.matrix_data.toarray(),
                               [[500, 200, 0], [500, 800, 0], [0, 0, 0]])
        self.assertEqual(list(out.ids()), ['s1', 's2', 'empty'])
        self.assertEqual(list(out.ids(axis='observation')), ['f1', 'f2', 'f3'])

    def test_total_sets_the_floor(self):
        # 4e-7 of s2 is lost at 1e6 (rounds to 0) and kept at 1e7
        lost = frequency(self.table, 1_000_000).get_value_by_ids('f3', 's2')
        kept = frequency(self.table, 10_000_000).get_value_by_ids('f3', 's2')
        self.assertEqual((lost, kept), (0, 4))

    def test_input_unchanged(self):
        before = self.table.copy()
        frequency(self.table, 1000)
        self.assertEqual(self.table, before)

    def test_not_relative(self):
        counts = biom.Table(np.array([[5.0, 2.0], [5.0, 8.0]]),
                            ['f1', 'f2'], ['s1', 's2'])
        with self.assertRaisesRegex(ValueError, 'do not sum to 1'):
            frequency(counts, 1000)

    def test_bad_total(self):
        with self.assertRaisesRegex(ValueError, 'at least 1'):
            frequency(self.table, 0)
