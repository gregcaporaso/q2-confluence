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
import numpy as np
import pandas as pd
import qiime2
import q2templates
from q2templates.reports import matryoshka_template

_ASSETS = importlib.resources.files('la_confluence') / 'assets'
EPS = 1e-9


def paired_report(ctx, table, tree, metadata, n_boots=10, permutations=999,
                  cv=10, n_estimators=500, random_state=None, n_jobs=1):
    """How far apart are the sources, on the samples every source profiled?

    `table` and `metadata` come from integrate (metadata columns
    original-sample-id and source). Samples whose original ID is present in
    every source are kept; the rest are dropped.
    """
    md = metadata.view(qiime2.Metadata).to_dataframe()
    for col in ('original-sample-id', 'source'):
        if col not in md.columns:
            raise ValueError(f'metadata has no {col!r} column; use the metadata '
                             'written by integrate')
    sources = sorted(md['source'].unique())
    by_source = [set(md.loc[md['source'] == s, 'original-sample-id'])
                 for s in sources]
    paired = set.intersection(*by_source)
    if len(paired) < 2:
        raise ValueError(f'{len(paired)} original sample IDs are present in every '
                         f'source ({sources}); at least 2 are needed')
    kept = md[md['original-sample-id'].isin(paired)]
    kept_md = qiime2.Metadata(kept)
    source_col = kept_md.get_column('source')

    filter_samples = ctx.get_action('feature_table', 'filter_samples')
    filter_features = ctx.get_action('feature_table', 'filter_features')
    core_metrics = ctx.get_action('boots', 'core_metrics')
    permanova = ctx.get_action('diversity', 'beta_group_significance')
    classify = ctx.get_action('sample_classifier', 'classify_samples_ncv')
    confusion = ctx.get_action('sample_classifier', 'confusion_matrix')
    heatmap = ctx.get_action('sample_classifier', 'heatmap')
    relative_frequency = ctx.get_action('feature_table', 'relative_frequency')
    summarize = ctx.get_action('confluence', 'summarize_paired_features')

    table, = filter_samples(table=table, metadata=kept_md)
    table, = filter_features(table=table, min_frequency=1)
    depth = int(table.view(biom.Table).sum(axis='sample').min())
    boots = core_metrics(table=table, phylogeny=tree, sampling_depth=depth,
                         metadata=kept_md, n=n_boots, replacement=True,
                         random_seed=random_state)
    sections = {'Ordination of the bootstrapped distances': boots.scatter_plot}
    for metric, dm in boots.distance_matrices.items():
        viz, = permanova(distance_matrix=dm, metadata=source_col,
                         method='permanova', pairwise=True,
                         permutations=permutations)
        sections[f'PERMANOVA by source: {metric.replace("_", " ")}'] = viz

    rf = classify(table=table, metadata=source_col, cv=cv,
                  random_state=random_state, n_jobs=n_jobs,
                  n_estimators=n_estimators, estimator='RandomForestClassifier')
    sections['Random forest: confusion matrix (nested cross-validation)'], = \
        confusion(predictions=rf.predictions, probabilities=rf.probabilities,
                  truth=source_col)
    sections['Random forest: heatmap of the most important features'] = \
        heatmap(table=table, importance=rf.feature_importance,
                sample_metadata=source_col).heatmap
    relative, = relative_frequency(table=table)
    sections['Per-feature summary'], = summarize(
        table=relative, importance=rf.feature_importance,
        predictions=rf.predictions, metadata=kept_md)
    report = ctx.make_report(matryoshka_template, sections)
    return report, boots.distance_matrices, rf.feature_importance


def summarize_paired_features(output_dir: str, table: biom.Table,
                              importance: pd.DataFrame, predictions: pd.Series,
                              metadata: qiime2.Metadata) -> None:
    """Importance next to each source's abundance and prevalence per feature."""
    md = metadata.to_dataframe()
    source = md['source']
    sources = sorted(source.unique())
    df = table.norm(axis='sample', inplace=False).to_dataframe(dense=True)
    out = pd.DataFrame(index=df.index)
    out.index.name = 'feature-id'
    out['importance'] = importance.reindex(df.index).iloc[:, 0].fillna(0.0)
    out['rank'] = out['importance'].rank(ascending=False, method='min').astype(int)
    for s in sources:
        cols = [c for c in df.columns if source.get(c) == s]
        out[f'mean-{s}'] = df[cols].mean(axis=1)
        out[f'prevalence-{s}'] = (df[cols] > 0).sum(axis=1).astype(int)
    a, b = sources[0], sources[1]
    out[f'log2-{b}-over-{a}'] = np.log2((out[f'mean-{b}'] + EPS)
                                        / (out[f'mean-{a}'] + EPS))
    only = {s: (out[[f'prevalence-{t}' for t in sources if t != s]] == 0).all(axis=1)
            & (out[f'prevalence-{s}'] > 0) for s in sources}
    out['direction'] = 'shared'
    for s in sources:
        out.loc[only[s], 'direction'] = f'{s}-only'
    out = out.sort_values('rank')
    out.to_csv(os.path.join(output_dir, 'feature-summary.tsv'), sep='\t')

    truth = source.reindex(predictions.index)
    accuracy = float((predictions == truth).mean())
    share = (out.groupby('direction')['importance'].sum()
             / out['importance'].sum()).rename('importance share').to_frame()
    share['features'] = out.groupby('direction').size()
    html = dict(classes='table table-striped', border=0, float_format='%.4g')
    q2templates.render(str(_ASSETS / 'paired_features' / 'index.html'), output_dir,
                       context={'accuracy': f'{accuracy:.2f}',
                                'n_samples': len(predictions),
                                'sources': ', '.join(sources),
                                'share': share.to_html(**html),
                                'top': out.head(50).to_html(**html)})
