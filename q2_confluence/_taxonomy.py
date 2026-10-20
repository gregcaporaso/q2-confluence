# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re

# A rank token with no name: '', '__', or a bare prefix such as 'g__'.
_UNNAMED = re.compile(r'^([a-z]__)?(__)?$')


def _collapsed_name(taxon: str, level: int, n_ranks: int) -> str:
    """The feature ID q2-taxa's `collapse` gives a taxon at `level`.

    Tokens are stripped, the string is padded with '__' to `n_ranks`
    tokens (the longest string in the taxonomy), and the first `level`
    tokens are joined with ';'. Keeping this rule identical to q2-taxa's
    is what lets collapsed tree tips match collapsed table features.

    The rule is the inner `_collapse` of `q2_taxa._util._collapse_table`
    (q2-taxa 2026.7.0), with `n_ranks` playing `max_observed_level`, which
    `q2_taxa._util._get_max_level` takes from the taxonomy passed to
    `collapse`; `test_taxonomy` runs that function on the same strings. The
    names agree as long as both taxonomies have at least `level` ranks,
    which holds for the seven-rank GTDB strings used throughout.
    """
    tokens = [t.strip() for t in taxon.split(';')]
    tokens += ['__'] * (n_ranks - len(tokens))
    return ';'.join(tokens[:level])


def _is_named(name: str) -> bool:
    """True when the last rank token of a collapsed name carries a name."""
    return _UNNAMED.match(name.rsplit(';', 1)[-1].strip()) is None
