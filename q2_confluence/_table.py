# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import numpy as np


def frequency(table: biom.Table, total: int) -> biom.Table:
    """Turn proportions into counts: round(proportion * total) per sample.

    Every sample must sum to 1 (empty samples, summing to 0, pass through).
    A sample's counts sum to `total` up to rounding, and a feature below
    0.5 / total in a sample rounds to zero there, so `total` sets the
    smallest proportion that survives.
    """
    if total < 1:
        raise ValueError(f'total must be at least 1, got {total}')
    sums = table.sum(axis='sample')
    bad = [s for s, v in zip(table.ids(), sums)
           if v != 0 and not np.isclose(v, 1.0, atol=1e-6)]
    if bad:
        raise ValueError(f'{len(bad)} samples do not sum to 1, e.g. {bad[:3]}; '
                         'the input must be a relative frequency table')
    return table.transform(lambda v, i, m: np.round(v * total),
                           axis='sample', inplace=False)
