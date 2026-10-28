# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib.resources
import os

import biom
import pandas as pd
import qiime2
import q2templates
import skbio

_ASSETS = importlib.resources.files('q2_confluence') / 'assets'


def summarize_inputs(output_dir: str, tables: biom.Table, tree: skbio.TreeNode,
                     metadata: qiime2.Metadata = None) -> None:
    # TODO: this is AI-generated code; human review pending
    """What integrate will find in its inputs: per study, the samples of
    its table and metadata and the features that are tips of the tree;
    across studies, the sample IDs the tables share (the pairs) and the
    feature IDs the tables and the tree share. `tables` and `metadata`
    are collections keyed by study, as for integrate."""
    metadata = metadata or {}
    unknown = sorted(set(metadata) - set(tables))
    if unknown:
        raise ValueError(f'metadata keys {unknown} are not table keys; '
                         'key each study\'s metadata like its table')
    tips = {t.name for t in tree.tips()}
    samples = {k: set(t.ids()) for k, t in tables.items()}
    features = {k: set(t.ids(axis='observation')) for k, t in tables.items()}
    rows = []
    for k in tables:
        ids = set(metadata[k].ids) if k in metadata else None
        in_tree = len(features[k] & tips)
        rows.append({
            'study': k, 'samples': len(samples[k]),
            'metadata-rows': len(ids) if ids is not None else '',
            'samples-with-metadata': len(samples[k] & ids) if ids is not None else '',
            'features': len(features[k]), 'features-in-tree': in_tree,
            'features-kept-share': round(in_tree / len(features[k]), 4)})
    studies = pd.DataFrame(rows)

    sample_ids = {f'table-{k}': samples[k] for k in tables}
    sample_ids |= {f'metadata-{k}': set(metadata[k].ids) for k in metadata}
    all_ids = sorted(set.union(*sample_ids.values()))
    membership = pd.DataFrame({k: [i in ids for i in all_ids]
                               for k, ids in sample_ids.items()},
                              index=pd.Index(all_ids, name='sample-id'))
    membership.to_csv(os.path.join(output_dir, 'sample-membership.tsv'), sep='\t')

    html = dict(classes='table table-striped', border=0)
    q2templates.render(str(_ASSETS / 'summarize_inputs' / 'index.html'), output_dir,
                       context={'studies': studies.to_html(index=False, **html),
                                'sample_overlap': _overlap(samples).to_html(**html),
                                'n_samples': len(all_ids),
                                'paired': len(set.intersection(*samples.values())),
                                'feature_overlap': _overlap(
                                    features | {'tree': tips}).to_html(**html),
                                'n_features': len(set.union(*features.values())),
                                'kept_everywhere': len(
                                    set.intersection(*features.values()) & tips)})


def _overlap(sets):
    """Square table: the number of IDs each pair of sets shares."""
    labels = list(sets)
    return pd.DataFrame([[len(sets[a] & sets[b]) for b in labels] for a in labels],
                        index=labels, columns=labels)
