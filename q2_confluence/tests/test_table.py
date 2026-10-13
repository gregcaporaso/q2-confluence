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

from q2_confluence._table import relative_frequency


class RelativeFrequencyTests(unittest.TestCase):

    def test_counts(self):
        # TODO: this is AI-generated code; human review pending
        table = biom.Table(np.array([[12.0, 0.0], [3.0, 9.0]]),
                           ['f1', 'f2'], ['s1', 's2'])
        out = relative_frequency(table)
        npt.assert_allclose(out.matrix_data.toarray(), [[0.8, 0.0], [0.2, 1.0]])
        self.assertEqual(list(out.ids()), ['s1', 's2'])
        self.assertEqual(list(out.ids(axis='observation')), ['f1', 'f2'])

    def test_proportions_after_a_drop(self):
        # TODO: this is AI-generated code; human review pending
        table = biom.Table(np.array([[0.6, 0.1], [0.3, 0.85]]),
                           ['f1', 'f2'], ['s1', 's2'])
        out = relative_frequency(table)
        npt.assert_allclose(out.matrix_data.toarray(),
                            [[2 / 3, 0.1 / 0.95], [1 / 3, 0.85 / 0.95]])

    def test_empty_sample(self):
        # TODO: this is AI-generated code; human review pending
        table = biom.Table(np.array([[1.0, 0.0], [1.0, 0.0]]),
                           ['f1', 'f2'], ['s1', 'empty'])
        out = relative_frequency(table)
        npt.assert_allclose(out.matrix_data.toarray(), [[0.5, 0.0], [0.5, 0.0]])

    def test_input_unchanged(self):
        # TODO: this is AI-generated code; human review pending
        table = biom.Table(np.array([[12.0, 0.0], [3.0, 9.0]]),
                           ['f1', 'f2'], ['s1', 's2'])
        before = table.copy()
        relative_frequency(table)
        self.assertEqual(table, before)
