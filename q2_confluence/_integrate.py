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


def integrate(ctx, tables, taxonomies, tree, metadata=None):
    # TODO: this is AI-generated code; human review pending
    """Merge every table in one feature space: the tips of `tree`.

    `tables` (counts or proportions, already keyed by the tree's tips, e.g.
    GTDB genome accessions) and `taxonomies` are collections with the same
    keys. Every source is renormalized to proportions and merged; the
    taxonomy of the merged features comes from the input taxonomies. Sample
    IDs in the merged table are '<id>-<key>', and the output metadata
    records the original ID (`original-sample-id`) and the key
    (`profile-type`), joined to the union of the studies' metadata when
    given.
    """
    # REVIEW COMMENT: the need to iterate explicitly over keys may have been
    # fixed in 2026.10. We'll install a 2026.10 environment and test doing this
    # without .keys() (I think that's a workaround you had from an earlier
    # issue we ran into with this).
    keys = list(tables.keys())   # a ResultCollection: iterate explicitly
    if set(keys) != set(taxonomies.keys()):
        raise ValueError('tables and taxonomies must have the same keys: '
                         f'{sorted(keys)} vs {sorted(taxonomies.keys())}')
    relative_frequency = ctx.get_action('confluence', 'relative_frequency')
    rename_ids = ctx.get_action('feature_table', 'rename_ids')
    merge = ctx.get_action('feature_table', 'merge')
    merge_taxa = ctx.get_action('feature_table', 'merge_taxa')
    # REVIEW COMMENT: why do we need our own merge metadata?
    merge_metadata = ctx.get_action('confluence', 'merge_metadata')

    relatives, rows = [], []
    for key in keys:
        relative, = relative_frequency(table=tables[key])
        new_ids = {i: f'{i}-{key}' for i in relative.view(biom.Table).ids()}
        renamed, = rename_ids(table=relative, metadata=_column(new_ids),
                              axis='sample', strict=True)
        relatives.append(renamed)
        rows += [{'sample-id': new, 'original-sample-id': old, 'profile-type': key}
                 for old, new in new_ids.items()]

    merged, = merge(tables=relatives)
    _check_tips(merged, tree)
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
    # Read like gtdb_tree_collapse, with underscores kept: the table's IDs
    # are the tips' raw names, as UniFrac's own parser also reads them.
    t = skbio.TreeNode.read(str(newick), format='newick',
                            convert_underscores=False)
    tips = {n.name for n in t.tips()}
    missing = [str(i) for i in table.view(biom.Table).ids(axis='observation')
               if i not in tips]
    if missing:
        raise ValueError(f'{len(missing)} features are not tips of the tree, e.g. '
                         f'{missing[:3]}; the tables must be keyed like the tree '
                         '(GTDB genome accessions, RS_GCF_... and GB_GCA_...), '
                         'from the same reference release')
