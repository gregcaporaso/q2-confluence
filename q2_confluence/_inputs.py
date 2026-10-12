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

_ASSETS = importlib.resources.files('q2_confluence') / 'assets'


def summarize_inputs(output_dir: str, tables: biom.Table,
                     metadata: qiime2.Metadata = None) -> None:
    """Which sample IDs do the inputs share? `tables` and `metadata` are
    collections (dicts keyed by label); the labels must not collide."""
    ids = {f'table-{k}': set(t.ids()) for k, t in tables.items()}
    for k, md in (metadata or {}).items():
        key = f'metadata-{k}'
        if key in ids:
            raise ValueError(f'duplicate input label {k!r}')
        ids[key] = set(md.ids)
    labels = list(ids)

    inputs = pd.DataFrame({
        'input': labels,
        'samples': [len(ids[k]) for k in labels],
        'features': [len(tables[k[6:]].ids(axis='observation'))
                     if k.startswith('table-') else '' for k in labels]})
    overlap = pd.DataFrame([[len(ids[a] & ids[b]) for b in labels] for a in labels],
                           index=labels, columns=labels)
    all_ids = sorted(set.union(*ids.values()))
    membership = pd.DataFrame({k: [i in ids[k] for i in all_ids] for k in labels},
                              index=pd.Index(all_ids, name='sample-id'))
    membership.to_csv(os.path.join(output_dir, 'sample-membership.tsv'), sep='\t')
    patterns = (membership.apply(lambda r: ', '.join(k for k in labels if r[k]), axis=1)
                .value_counts().rename_axis('present-in').reset_index(name='samples'))
    in_all = sum(membership.all(axis=1))

    html = dict(index=False, classes='table table-striped', border=0)
    q2templates.render(str(_ASSETS / 'summarize_inputs' / 'index.html'), output_dir,
                       context={'inputs': inputs.to_html(**html),
                                'overlap': overlap.to_html(
                                    classes='table table-striped', border=0),
                                'patterns': patterns.to_html(**html),
                                'n_all': len(all_ids), 'in_all': in_all})
