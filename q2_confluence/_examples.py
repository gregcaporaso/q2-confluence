# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Usage examples, registered on the actions and rendered by q2doc.

The factories build small artifacts shaped like the real inputs: a GTDB-style
tree with bootstrap:taxon labels on internal nodes and genome accessions as
tips, and the matching genome-keyed taxonomy.
"""


import pandas as pd
from qiime2 import Artifact, Metadata
from q2_confluence._formats import GTDBNewickFormat

_TREE = ("(((RS_GCF_1:0.10,RS_GCF_2:0.30)'100.0:g__Blautia':0.50,"
         "GB_GCA_3:0.40)'90.0:f__Lachnospiraceae':1.00,"
         "(RS_GCF_4:0.20,RS_GCF_5:0.20)'95.0:g__Bacillus':2.00);")

_TAXONOMY = {
    'RS_GCF_1': 'd__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
                'f__Lachnospiraceae;g__Blautia;s__Blautia wexlerae',
    'RS_GCF_2': 'd__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
                'f__Lachnospiraceae;g__Blautia;s__Blautia_A obeum',
    'GB_GCA_3': 'd__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
                'f__Lachnospiraceae;g__Roseburia;s__Roseburia hominis',
    'RS_GCF_4': 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;'
                'f__Bacillaceae;g__Bacillus;s__Bacillus subtilis',
    'RS_GCF_5': 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;'
                'f__Bacillaceae;g__Bacillus;s__Bacillus licheniformis',
}


def _reference_tree_factory():
    # TODO: this is AI-generated code; human review pending
    # Imported as a GTDB tree is: the tips get quoted on the way in.
    fmt = GTDBNewickFormat()
    with fmt.open() as fh:
        fh.write(_TREE)
    return Artifact.import_data('Phylogeny[Rooted]', fmt)


def _reference_taxonomy_factory():
    # TODO: this is AI-generated code; human review pending
    df = pd.DataFrame({'Taxon': pd.Series(_TAXONOMY)})
    df.index.name = 'Feature ID'
    return Artifact.import_data('FeatureData[Taxonomy]', df)


def _gtdb_tree_collapse_genus(use):
    # TODO: this is AI-generated code; human review pending
    tree = use.init_artifact('reference-tree', _reference_tree_factory)
    taxonomy = use.init_artifact('reference-taxonomy',
                                 _reference_taxonomy_factory)
    use.action(
        use.UsageAction('confluence', 'gtdb_tree_collapse'),
        use.UsageInputs(tree=tree, taxonomy=taxonomy, level=6),
        use.UsageOutputNames(collapsed_tree='genus-tree',
                             collapse_log='genus-tree-collapse-log'),
    )


def _relative_table_factory():
    # TODO: this is AI-generated code; human review pending
    import biom
    import numpy as np
    # proportions left by dropping a feature: the samples sum to 0.9 and 0.95
    table = biom.Table(np.array([[0.6, 0.1], [0.3, 0.85]]),
                       ['d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
                        'f__Lachnospiraceae;g__Blautia',
                        'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;'
                        'f__Bacillaceae;g__Bacillus'],
                       ['sample-1', 'sample-2'])
    return Artifact.import_data('FeatureTable[RelativeFrequency]', table)


def _relative_frequency_example(use):
    # TODO: this is AI-generated code; human review pending
    table = use.init_artifact('filtered-relative-table', _relative_table_factory)
    use.action(
        use.UsageAction('confluence', 'relative_frequency'),
        use.UsageInputs(table=table),
        use.UsageOutputNames(relative_frequency_table='relative-table'),
    )


def _amplicon_table_factory():
    # TODO: this is AI-generated code; human review pending
    import biom
    import numpy as np
    table = biom.Table(np.array([[12, 0, 7], [3, 9, 1]]), ['asv-1', 'asv-2'],
                       ['sample-1', 'sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[Frequency]', table)


def _metagenome_table_factory():
    # TODO: this is AI-generated code; human review pending
    import biom
    import numpy as np
    table = biom.Table(np.array([[0.6, 0.1], [0.4, 0.9]]), ['GCF_1', 'GCF_2'],
                       ['sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[RelativeFrequency]', table)


def _metadata_factory(ids, **columns):
    # TODO: this is AI-generated code; human review pending
    df = pd.DataFrame(columns, index=pd.Index(ids, name='sample-id'))
    return Artifact.import_data('ImmutableMetadata', Metadata(df))


def _amplicon_metadata_factory():
    # TODO: this is AI-generated code; human review pending
    return _metadata_factory(['sample-1', 'sample-2', 'sample-3'],
                             site=['compost', 'compost', 'soil'],
                             ph=[6.8, 7.2, 7.9])


def _metagenome_metadata_factory():
    # TODO: this is AI-generated code; human review pending
    return _metadata_factory(['sample-2', 'sample-3', 'sample-4'],
                             site=['compost', 'soil', 'soil'],
                             **{'read-pairs-millions': [41.2, 38.7, 27.5]})


def _summarize_inputs_two_studies(use):
    # TODO: this is AI-generated code; human review pending
    amplicon = use.init_artifact('amplicon-table', _amplicon_table_factory)
    metagenome = use.init_artifact('metagenome-table', _metagenome_table_factory)
    amplicon_md = use.init_artifact('amplicon-metadata', _amplicon_metadata_factory)
    metagenome_md = use.init_artifact('metagenome-metadata',
                                      _metagenome_metadata_factory)
    tables = use.construct_artifact_collection(
        'tables', {'amplicon': amplicon, 'metagenome': metagenome})
    metadata = use.construct_artifact_collection(
        'metadata', {'amplicon': amplicon_md, 'metagenome': metagenome_md})
    use.action(
        use.UsageAction('confluence', 'summarize_inputs'),
        use.UsageInputs(tables=tables, metadata=metadata),
        use.UsageOutputNames(visualization='inputs-summary'),
    )


def _taxonomy_factory(mapping):
    # TODO: this is AI-generated code; human review pending
    df = pd.DataFrame({'Taxon': pd.Series(mapping)})
    df.index.name = 'Feature ID'
    return Artifact.import_data('FeatureData[Taxonomy]', df)


# --- the genome space: tables keyed by the reference tree's tips -------------

_GENOMES = ['RS_GCF_1', 'GB_GCA_3', 'RS_GCF_4']


def _amplicon_genome_table_factory():
    # TODO: this is AI-generated code; human review pending
    """ASV counts summed per reference genome, as closed-reference
    clustering against the GTDB SSU reference leaves them."""
    import biom
    import numpy as np
    table = biom.Table(np.array([[12, 0, 7], [3, 9, 1], [1, 2, 0]]), _GENOMES,
                       ['sample-1', 'sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[Frequency]', table)


def _metagenome_genome_table_factory():
    # TODO: this is AI-generated code; human review pending
    """sylph's taxonomic abundance per genome, GTDB-prefixed accessions."""
    import biom
    import numpy as np
    table = biom.Table(np.array([[0.6, 0.1], [0.4, 0.9]]),
                       ['RS_GCF_1', 'RS_GCF_2'], ['sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[RelativeFrequency]', table)


def _integrate_genomes(use):
    # TODO: this is AI-generated code; human review pending
    tree = use.init_artifact('reference-tree', _reference_tree_factory)
    amplicon = use.init_artifact('amplicon-genome-table',
                                 _amplicon_genome_table_factory)
    metagenome = use.init_artifact('metagenome-genome-table',
                                   _metagenome_genome_table_factory)
    amplicon_md = use.init_artifact('amplicon-metadata', _amplicon_metadata_factory)
    metagenome_md = use.init_artifact('metagenome-metadata',
                                      _metagenome_metadata_factory)
    tables = use.construct_artifact_collection(
        'tables', {'amplicon': amplicon, 'metagenome': metagenome})
    metadata = use.construct_artifact_collection(
        'metadata', {'amplicon': amplicon_md, 'metagenome': metagenome_md})
    use.action(
        use.UsageAction('confluence', 'integrate'),
        use.UsageInputs(tables=tables, tree=tree, metadata=metadata),
        use.UsageOutputNames(integrated_table='genome-relative-table',
                             integrated_metadata='genome-metadata'),
    )


# --- paired study: six samples profiled by both sources ----------------------

_PAIRED_SAMPLES = [f'sample-{i}' for i in range(1, 7)]


def _paired_counts(rng, bias):
    # TODO: this is AI-generated code; human review pending
    """Six samples x three taxa; `bias` tilts the third taxon by source so the
    random forest has a signal."""
    base = rng.integers(50, 500, size=(3, 6)).astype(float)
    base[2] *= bias
    return base


def _paired_amplicon_table_factory():
    # TODO: this is AI-generated code; human review pending
    """Closed-reference counts per genome; the same three genomes sylph sees."""
    import biom
    import numpy as np
    rng = np.random.default_rng(1)
    table = biom.Table(_paired_counts(rng, 0.2), _GENOMES, _PAIRED_SAMPLES)
    return Artifact.import_data('FeatureTable[Frequency]', table)


def _paired_metagenome_table_factory():
    # TODO: this is AI-generated code; human review pending
    import biom
    import numpy as np
    rng = np.random.default_rng(2)
    counts = _paired_counts(rng, 5.0) * 1000
    table = biom.Table(counts, _GENOMES, _PAIRED_SAMPLES)
    return Artifact.import_data('FeatureTable[Frequency]', table)


def _paired_metadata_factory():
    # TODO: this is AI-generated code; human review pending
    return _metadata_factory(_PAIRED_SAMPLES,
                             site=['compost'] * 3 + ['soil'] * 3)


def _paired_report_example(use):
    # TODO: this is AI-generated code; human review pending
    tree = use.init_artifact('reference-tree', _reference_tree_factory)
    taxonomy = use.init_artifact('reference-taxonomy',
                                 _reference_taxonomy_factory)
    tables = use.construct_artifact_collection('tables', {
        'amplicon': use.init_artifact('amplicon-genome-table',
                                      _paired_amplicon_table_factory),
        'metagenome': use.init_artifact('metagenome-genome-table',
                                        _paired_metagenome_table_factory)})
    study = use.init_artifact('study-metadata', _paired_metadata_factory)
    # one study profiled twice: its metadata goes in under both keys
    metadata = use.construct_artifact_collection(
        'metadata', {'amplicon': study, 'metagenome': study})
    integrated_table, merged_metadata = use.action(
        use.UsageAction('confluence', 'integrate'),
        use.UsageInputs(tables=tables, tree=tree, metadata=metadata),
        use.UsageOutputNames(integrated_table='genome-relative-table',
                             integrated_metadata='genome-metadata'),
    )
    use.action(
        use.UsageAction('confluence', 'paired_report'),
        use.UsageInputs(table=integrated_table, tree=tree,
                        metadata=merged_metadata, taxonomy=taxonomy,
                        permutations=99, cv=3, n_estimators=50,
                        random_state=1),
        use.UsageOutputNames(report='paired-report',
                             distance_matrices='distance-matrices',
                             feature_importance='feature-importance',
                             feature_summary='feature-summary'),
    )


# --- fixtures of tests/test_report.py for the per-feature summary: a merged
# --- relative table, its metadata, and the random forest's importances, as
# --- integrate and paired-report would produce them for three features of
# --- two studies (the hidden action has no usage example)

_PAIRED_FEATURES = _GENOMES
_PAIRED_MERGED_IDS = [f'{s}-{k}' for k in ('amplicon', 'metagenome')
                      for s in _PAIRED_SAMPLES[:3]]


def _paired_relative_table_factory():
    # TODO: this is AI-generated code; human review pending
    import biom
    import numpy as np
    data = np.array([[0.5, 0.6, 0.4, 0.3, 0.2, 0.3],    # shared
                     [0.5, 0.4, 0.6, 0.0, 0.0, 0.0],    # amplicon only
                     [0.0, 0.0, 0.0, 0.7, 0.8, 0.7]])   # metagenome only
    table = biom.Table(data, _PAIRED_FEATURES, _PAIRED_MERGED_IDS)
    return Artifact.import_data('FeatureTable[RelativeFrequency]', table)


def _paired_integrated_metadata_factory():
    # TODO: this is AI-generated code; human review pending
    df = pd.DataFrame({'original-sample-id': _PAIRED_SAMPLES[:3] * 2,
                       'study-id': ['amplicon'] * 3 + ['metagenome'] * 3},
                      index=pd.Index(_PAIRED_MERGED_IDS, name='sample-id'))
    return Artifact.import_data('ImmutableMetadata', Metadata(df))


def _paired_importance_factory():
    # TODO: this is AI-generated code; human review pending
    df = pd.DataFrame({'importance': [0.1, 0.5, 0.4]},
                      index=pd.Index(_PAIRED_FEATURES, name='feature-id'))
    return Artifact.import_data('FeatureData[Importance]', df)
