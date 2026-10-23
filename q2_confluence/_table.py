# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def relative_frequency(table: biom.Table) -> biom.Table:
    """Divide every sample by its total, as feature-table relative-frequency
    does, but from counts or from proportions: a relative table from which
    features were dropped sums to 1 again. An empty sample stays empty."""
    return table.norm(axis='sample', inplace=False)
