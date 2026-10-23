# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import warnings

import biom
import pandas as pd
import qiime2
from rachis.core.exceptions import RachisWarning


# A source that loses more than this share of its features to the tree
# filter gets a visible warning (its counts are in the metadata either way).
DROPPED_SHARE_WARNING = 0.1


def integrate(ctx, tables, tree, metadata=None):
    # TODO: this is AI-generated code; human review pending (changed after
    # review on 2026-10-23: study-id column, metadata keys checked first,
    # collections iterated directly)
    """Merge every table in one feature space: the tips of `tree`.

    `tables` (counts or proportions, keyed by the tree's tips, e.g. GTDB
    genome accessions) is a collection keyed by a label. For each source the
    features that are not tips of the tree are dropped (phylogeny
    filter-table) and counted, the rest is renormalized to proportions and
    the sample IDs are suffixed with the key; the sources are then merged.
    The output metadata records, per merged sample, `original-sample-id`,
    `study-id` (the key), the source's `features-before` and
    `features-after` the tree filter and their quotient
    `features-kept-share`, and the columns of the sample's own study
    metadata when given: `metadata` is keyed like `tables`, each study's
    sample IDs are suffixed with its key, and the studies are then merged
    with metadata merge, which they pass because no two share an ID. A
    source that loses more than DROPPED_SHARE_WARNING of its features
    raises a visible warning; one that loses every feature is an error.
    """
    if metadata:
        # a usage error should fail before any work is done
        unknown = sorted(set(metadata) - set(tables))
        if unknown:
            raise ValueError(f'metadata keys {unknown} are not table keys; '
                             'key each study\'s metadata like its table')
    filter_table = ctx.get_action('phylogeny', 'filter_table')
    relative_frequency = ctx.get_action('confluence', 'relative_frequency')
    rename_ids = ctx.get_action('feature_table', 'rename_ids')
    merge = ctx.get_action('feature_table', 'merge')
    merge_metadata = ctx.get_action('metadata', 'merge')

    relatives, rows = [], []
    for key in tables:
        before = len(tables[key].view(biom.Table).ids(axis='observation'))
        kept, = filter_table(table=tables[key], tree=tree)
        after = len(kept.view(biom.Table).ids(axis='observation'))
        if after == 0:
            raise ValueError(
                f'{key}: none of its {before} features is a tip of the tree; '
                'the tables must be keyed like the tree (GTDB genome '
                'accessions, RS_GCF_... and GB_GCA_...), from the same '
                'reference release')
        share = after / before
        if 1 - share > DROPPED_SHARE_WARNING:
            warnings.warn(
                f'{key}: {before - after} of its {before} features '
                f'({1 - share:.0%}) are not tips of the tree and were '
                'dropped. Genomes outside the tree (archaea against the '
                'bacterial tree) are expected; a large share suggests a '
                'reference release mismatch. The counts are in the '
                'integrated metadata (features-before, features-after, '
                'features-kept-share).', RachisWarning)
        relative, = relative_frequency(table=kept)
        new_ids = {i: f'{i}-{key}' for i in relative.view(biom.Table).ids()}
        renamed, = rename_ids(table=relative, metadata=_column(new_ids),
                              axis='sample', strict=True)
        relatives.append(renamed)
        rows += [{'sample-id': new, 'original-sample-id': old,
                  'study-id': key, 'features-before': before,
                  'features-after': after, 'features-kept-share': share}
                 for old, new in new_ids.items()]

    merged, = merge(tables=relatives)

    md = pd.DataFrame(rows).set_index('sample-id')
    if metadata:
        studies = [_suffixed(metadata[key].view(qiime2.Metadata), key)
                   for key in metadata]
        union = studies[0]
        for other in studies[1:]:
            merged_md, = merge_metadata(metadata1=union, metadata2=other,
                                        method='union')
            union = merged_md.view(qiime2.Metadata)
        md = md.join(union.to_dataframe())
    integrated = ctx.make_artifact('ImmutableMetadata', qiime2.Metadata(md))
    return merged, integrated


def _suffixed(study, key):
    """The study's metadata with every sample ID suffixed by `-<key>`."""
    df = study.to_dataframe()
    df.index = pd.Index([f'{i}-{key}' for i in df.index], name='sample-id')
    return qiime2.Metadata(df)


def _column(mapping):
    s = pd.Series(mapping, name='new-id')
    s.index.name = 'sample-id'
    return qiime2.CategoricalMetadataColumn(s)
