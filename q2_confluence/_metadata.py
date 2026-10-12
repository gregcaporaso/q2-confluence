# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2


def merge_metadata(metadata: qiime2.Metadata,
                   how: str = 'union') -> qiime2.Metadata:
    """Merge the metadata of several studies into one.

    `metadata` is a collection (dict keyed by label). With 'union' every ID
    in any input is kept; with 'intersection' only the IDs in every input.
    Columns with the same name are combined; an ID with two different
    non-missing values in such a column is an error.
    """
    frames = {label: md.to_dataframe() for label, md in metadata.items()}
    merged = _merge_frames(list(frames.values()))
    if how == 'intersection':
        shared = set.intersection(*(set(f.index) for f in frames.values()))
        if not shared:
            raise ValueError('no sample ID is present in every metadata input; '
                             'nothing to intersect')
        merged = merged.loc[sorted(shared)]
    merged.index.name = 'sample-id'
    return qiime2.Metadata(merged)


def _merge_frames(frames):
    """Outer join on the index; shared columns combined, conflicts rejected."""
    out = frames[0]
    for f in frames[1:]:
        for col in out.columns.intersection(f.columns):
            ids = out.index.intersection(f.index)
            a, b = out.loc[ids, col], f.loc[ids, col]
            clash = a.notna() & b.notna() & (a.astype(str) != b.astype(str))
            if clash.any():
                raise ValueError(f'column {col!r} has different values for '
                                 f'{int(clash.sum())} shared IDs, e.g. '
                                 f'{list(ids[clash][:3])}')
        out = out.combine_first(f)
    return out.sort_index()
