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
from q2_types.tree import NewickFormat

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


def reference_tree_factory():
    fmt = NewickFormat()
    with fmt.open() as fh:
        fh.write(_TREE)
    return Artifact.import_data('Phylogeny[Rooted]', fmt)


def reference_taxonomy_factory():
    df = pd.DataFrame({'Taxon': pd.Series(_TAXONOMY)})
    df.index.name = 'Feature ID'
    return Artifact.import_data('FeatureData[Taxonomy]', df)


def gtdb_tree_collapse_genus(use):
    tree = use.init_artifact('reference-tree', reference_tree_factory)
    taxonomy = use.init_artifact('reference-taxonomy',
                                 reference_taxonomy_factory)
    use.action(
        use.UsageAction('confluence', 'gtdb_tree_collapse'),
        use.UsageInputs(tree=tree, taxonomy=taxonomy, level=6),
        use.UsageOutputNames(collapsed_tree='genus-tree'),
    )


def relative_table_factory():
    import biom
    import numpy as np
    table = biom.Table(np.array([[0.6, 0.1], [0.4, 0.9]]),
                       ['d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
                        'f__Lachnospiraceae;g__Blautia',
                        'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;'
                        'f__Bacillaceae;g__Bacillus'],
                       ['sample-1', 'sample-2'])
    return Artifact.import_data('FeatureTable[RelativeFrequency]', table)


def frequency_one_million(use):
    table = use.init_artifact('relative-table', relative_table_factory)
    use.action(
        use.UsageAction('confluence', 'frequency'),
        use.UsageInputs(table=table, total=1_000_000),
        use.UsageOutputNames(table='frequency-table'),
    )


def amplicon_table_factory():
    import biom
    import numpy as np
    table = biom.Table(np.array([[12, 0, 7], [3, 9, 1]]), ['asv-1', 'asv-2'],
                       ['sample-1', 'sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[Frequency]', table)


def metagenome_table_factory():
    import biom
    import numpy as np
    table = biom.Table(np.array([[0.6, 0.1], [0.4, 0.9]]), ['GCF_1', 'GCF_2'],
                       ['sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[RelativeFrequency]', table)


def _metadata_factory(ids, **columns):
    df = pd.DataFrame(columns, index=pd.Index(ids, name='sample-id'))
    return Artifact.import_data('ImmutableMetadata', Metadata(df))


def amplicon_metadata_factory():
    return _metadata_factory(['sample-1', 'sample-2', 'sample-3'],
                             site=['compost', 'compost', 'soil'],
                             ph=[6.8, 7.2, 7.9])


def metagenome_metadata_factory():
    return _metadata_factory(['sample-2', 'sample-3', 'sample-4'],
                             site=['compost', 'soil', 'soil'],
                             read_pairs_millions=[41.2, 38.7, 27.5])


def summarize_inputs_two_studies(use):
    amplicon = use.init_artifact('amplicon-table', amplicon_table_factory)
    metagenome = use.init_artifact('metagenome-table', metagenome_table_factory)
    amplicon_md = use.init_artifact('amplicon-metadata', amplicon_metadata_factory)
    metagenome_md = use.init_artifact('metagenome-metadata',
                                      metagenome_metadata_factory)
    use.action(
        use.UsageAction('confluence', 'summarize_inputs'),
        use.UsageInputs(tables={'amplicon': amplicon, 'metagenome': metagenome},
                        metadata={'amplicon': amplicon_md,
                                  'metagenome': metagenome_md}),
        use.UsageOutputNames(visualization='inputs-summary'),
    )


def merge_metadata_union(use):
    amplicon_md = use.init_artifact('amplicon-metadata', amplicon_metadata_factory)
    metagenome_md = use.init_artifact('metagenome-metadata',
                                      metagenome_metadata_factory)
    use.action(
        use.UsageAction('confluence', 'merge_metadata'),
        use.UsageInputs(metadata={'amplicon': amplicon_md,
                                  'metagenome': metagenome_md},
                        how='union'),
        use.UsageOutputNames(merged_metadata='merged-metadata'),
    )
