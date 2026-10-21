# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
import qiime2
import skbio
from q2_types.tree import NewickFormat

from ._taxonomy import _is_named


def integrate(ctx, tables, taxonomies, tree, level=None, metadata=None):
    # TODO: this is AI-generated code; human review pending
    """Merge every table in one feature space: the tree's tips.

    `tables` (counts or proportions) and `taxonomies` are collections with
    the same keys. With `level`, every table is collapsed to it, what is
    unnamed there is dropped, and `tree` is the reference tree collapsed to
    the same level; without, the tables are taken as they are (already
    keyed by the tree's tips, e.g. genome accessions). Every source is
    renormalized to proportions. Sample IDs in the merged table are
    '<id>-<key>', and the output metadata records the original ID
    (`original-sample-id`) and the key (`profile-type`), joined to the union
    of the studies' metadata when given.
    """
    keys = list(tables.keys())   # a ResultCollection: iterate explicitly
    if set(keys) != set(taxonomies.keys()):
        raise ValueError('tables and taxonomies must have the same keys: '
                         f'{sorted(keys)} vs {sorted(taxonomies.keys())}')
    collapse = ctx.get_action('taxa', 'collapse')
    filter_features = ctx.get_action('feature_table', 'filter_features')
    relative_frequency = ctx.get_action('confluence', 'relative_frequency')
    rename_ids = ctx.get_action('feature_table', 'rename_ids')
    merge = ctx.get_action('feature_table', 'merge')
    merge_taxa = ctx.get_action('feature_table', 'merge_taxa')
    ids_to_taxonomy = ctx.get_action('taxa', 'feature_ids_to_taxonomy')
    merge_metadata = ctx.get_action('confluence', 'merge_metadata')

    relatives, rows = [], []
    for key in keys:
        kept = tables[key]
        if level is not None:
            kept, = collapse(table=kept, taxonomy=taxonomies[key], level=level)
            unnamed = [i for i in kept.view(biom.Table).ids(axis='observation')
                       if not _is_named(i)]
            if unnamed:
                kept, = filter_features(table=kept, metadata=_ids(unnamed),
                                        exclude_ids=True)
        relative, = relative_frequency(table=kept)
        new_ids = {i: f'{i}-{key}' for i in kept.view(biom.Table).ids()}
        renamed, = rename_ids(table=relative, metadata=_column(new_ids),
                              axis='sample', strict=True)
        relatives.append(renamed)
        rows += [{'sample-id': new, 'original-sample-id': old, 'profile-type': key}
                 for old, new in new_ids.items()]

    merged, = merge(tables=relatives)
    _check_tips(merged, tree)
    if level is not None:
        taxonomy, = ids_to_taxonomy(table=merged)
    else:
        all_taxa, = merge_taxa(data=[taxonomies[k] for k in keys])
        taxonomy = ctx.make_artifact(
            'FeatureData[Taxonomy]', _restrict_taxonomy(all_taxa, merged))

    md = pd.DataFrame(rows).set_index('sample-id')
    if metadata:
        studies, = merge_metadata(metadata=metadata, how='union')
        studies = studies.view(qiime2.Metadata).to_dataframe()
        md = md.join(studies, on='original-sample-id')
    merged_md = ctx.make_artifact('ImmutableMetadata', qiime2.Metadata(md))
    return merged, taxonomy, merged_md


def _ids(ids):
    # TODO: this is AI-generated code; human review pending
    return qiime2.Metadata(pd.DataFrame(index=pd.Index(ids, name='feature-id')))


def _column(mapping):
    # TODO: this is AI-generated code; human review pending
    s = pd.Series(mapping, name='new-id')
    s.index.name = 'sample-id'
    return qiime2.CategoricalMetadataColumn(s)


def _restrict_taxonomy(taxonomy, table):
    # TODO: this is AI-generated code; human review pending
    """The taxonomy's rows for the table's features, in table order."""
    features = list(table.view(biom.Table).ids(axis='observation'))
    taxa = taxonomy.view(pd.Series)
    missing = [str(i) for i in features if i not in taxa.index]
    if missing:
        raise ValueError(f'{len(missing)} features have no taxonomy, e.g. '
                         f'{missing[:3]}')
    out = taxa.reindex(features)
    out.index.name = 'Feature ID'
    return out


def _check_tips(table, tree):
    # TODO: this is AI-generated code; human review pending
    newick = tree.view(NewickFormat)
    t = skbio.TreeNode.read(str(newick), format='newick',
                            convert_underscores=False, verify=False)
    tips = {n.name for n in t.tips()}
    missing = [str(i) for i in table.view(biom.Table).ids(axis='observation')
               if i not in tips]
    if missing:
        raise ValueError(f'{len(missing)} features are not tips of the tree, e.g. '
                         f'{missing[:3]}; collapse the reference tree to the same '
                         'level, from the same reference release')
