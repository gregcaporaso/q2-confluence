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

import gzip
import os
import shutil
import tempfile

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
                             **{'read-pairs-millions': [41.2, 38.7, 27.5]})


def summarize_inputs_two_studies(use):
    amplicon = use.init_artifact('amplicon-table', amplicon_table_factory)
    metagenome = use.init_artifact('metagenome-table', metagenome_table_factory)
    amplicon_md = use.init_artifact('amplicon-metadata', amplicon_metadata_factory)
    metagenome_md = use.init_artifact('metagenome-metadata',
                                      metagenome_metadata_factory)
    tables = use.construct_artifact_collection(
        'tables', {'amplicon': amplicon, 'metagenome': metagenome})
    metadata = use.construct_artifact_collection(
        'metadata', {'amplicon': amplicon_md, 'metagenome': metagenome_md})
    use.action(
        use.UsageAction('confluence', 'summarize_inputs'),
        use.UsageInputs(tables=tables, metadata=metadata),
        use.UsageOutputNames(visualization='inputs-summary'),
    )


def merge_metadata_union(use):
    amplicon_md = use.init_artifact('amplicon-metadata', amplicon_metadata_factory)
    metagenome_md = use.init_artifact('metagenome-metadata',
                                      metagenome_metadata_factory)
    metadata = use.construct_artifact_collection(
        'metadata', {'amplicon': amplicon_md, 'metagenome': metagenome_md})
    use.action(
        use.UsageAction('confluence', 'merge_metadata'),
        use.UsageInputs(metadata=metadata, how='union'),
        use.UsageOutputNames(merged_metadata='merged-metadata'),
    )


_ASV_TAXONOMY = {
    'asv-1': 'd__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
             'f__Lachnospiraceae;g__Blautia;s__',
    'asv-2': 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;'
             'f__Bacillaceae;g__Bacillus;s__Bacillus subtilis',
    'asv-3': 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;'
             'f__Bacillaceae;g__;s__',
}

_GENOME_TAXONOMY = {
    'GCF_1': _TAXONOMY['RS_GCF_1'],
    'GCF_2': _TAXONOMY['RS_GCF_4'],
}


def _taxonomy_factory(mapping):
    df = pd.DataFrame({'Taxon': pd.Series(mapping)})
    df.index.name = 'Feature ID'
    return Artifact.import_data('FeatureData[Taxonomy]', df)


def amplicon_taxonomy_factory():
    return _taxonomy_factory(_ASV_TAXONOMY)


def metagenome_taxonomy_factory():
    return _taxonomy_factory(_GENOME_TAXONOMY)


def amplicon_table_three_asvs_factory():
    import biom
    import numpy as np
    table = biom.Table(np.array([[12, 0, 7], [3, 9, 1], [1, 2, 0]]),
                       ['asv-1', 'asv-2', 'asv-3'],
                       ['sample-1', 'sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[Frequency]', table)


def metagenome_counts_factory():
    import biom
    import numpy as np
    table = biom.Table(np.array([[600000, 100000], [400000, 900000]]),
                       ['GCF_1', 'GCF_2'], ['sample-2', 'sample-3'])
    return Artifact.import_data('FeatureTable[Frequency]', table)


def integrate_at_genus(use):
    tree = use.init_artifact('reference-tree', reference_tree_factory)
    taxonomy = use.init_artifact('reference-taxonomy',
                                 reference_taxonomy_factory)
    genus_tree, = use.action(
        use.UsageAction('confluence', 'gtdb_tree_collapse'),
        use.UsageInputs(tree=tree, taxonomy=taxonomy, level=6),
        use.UsageOutputNames(collapsed_tree='genus-tree'),
    )
    amplicon = use.init_artifact('amplicon-table',
                                 amplicon_table_three_asvs_factory)
    amplicon_tax = use.init_artifact('amplicon-taxonomy',
                                     amplicon_taxonomy_factory)
    metagenome = use.init_artifact('metagenome-table', metagenome_counts_factory)
    metagenome_tax = use.init_artifact('metagenome-taxonomy',
                                       metagenome_taxonomy_factory)
    amplicon_md = use.init_artifact('amplicon-metadata', amplicon_metadata_factory)
    metagenome_md = use.init_artifact('metagenome-metadata',
                                      metagenome_metadata_factory)
    tables = use.construct_artifact_collection(
        'tables', {'amplicon': amplicon, 'metagenome': metagenome})
    taxonomies = use.construct_artifact_collection(
        'taxonomies', {'amplicon': amplicon_tax, 'metagenome': metagenome_tax})
    metadata = use.construct_artifact_collection(
        'metadata', {'amplicon': amplicon_md, 'metagenome': metagenome_md})
    use.action(
        use.UsageAction('confluence', 'integrate'),
        use.UsageInputs(tables=tables, taxonomies=taxonomies, tree=genus_tree,
                        level=6, total=1_000_000, metadata=metadata),
        use.UsageOutputNames(relative_table='genus-relative-table',
                             frequency_table='genus-frequency-table',
                             taxonomy='genus-taxonomy',
                             metadata='genus-metadata'),
    )


# --- sylph: synthetic genomes, a database sketched from them, and reads -----

def _random_genome(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def _write_fasta(path, name, sequence):
    with open(path, 'w') as fh:
        fh.write(f'>{name}\n')
        for i in range(0, len(sequence), 80):
            fh.write(sequence[i:i + 80] + '\n')


def _write_reads(path_prefix, rng, genome, n, read_length=150, fragment=400):
    """Paired reads from random fragments of `genome` (R2 reverse complemented)."""
    comp = str.maketrans('ACGT', 'TGCA')
    with gzip.open(f'{path_prefix}_R1.fastq.gz', 'wt') as r1, \
            gzip.open(f'{path_prefix}_R2.fastq.gz', 'wt') as r2:
        for i in range(n):
            start = rng.randrange(0, len(genome) - fragment)
            frag = genome[start:start + fragment]
            r1.write(f'@read{i}/1\n{frag[:read_length]}\n+\n{"I" * read_length}\n')
            rev = frag[-read_length:][::-1].translate(comp)
            r2.write(f'@read{i}/2\n{rev}\n+\n{"I" * read_length}\n')


def synthetic_sylph_data(directory, seed=1):
    """Two 60 kb genomes, a sylph database of them with a two-genome
    taxonomy, and three paired-end samples: genome 1 only, both genomes, and
    reads from a third genome the database does not hold. Returns the
    database directory and the reads manifest path."""
    import random
    import subprocess

    rng = random.Random(seed)
    genomes = {'GCF_000001.1': _random_genome(rng, 60000),
               'GCF_000002.1': _random_genome(rng, 60000)}
    other = _random_genome(rng, 60000)
    gdir = os.path.join(directory, 'genomes')
    os.mkdir(gdir)
    for acc, seq in genomes.items():
        _write_fasta(os.path.join(gdir, f'{acc}_genomic.fna'), acc, seq)
    db = os.path.join(directory, 'db')
    os.mkdir(db)
    subprocess.run(['sylph', 'sketch', '-t', '1', '-g']
                   + [os.path.join(gdir, f) for f in sorted(os.listdir(gdir))]
                   + ['-o', os.path.join(db, 'database')],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    with gzip.open(os.path.join(db, 'taxonomy.tsv.gz'), 'wt') as fh:
        fh.write('GCF_000001.1\t' + _TAXONOMY['RS_GCF_4'] + '\n')
        fh.write('GCF_000002.1\t' + _TAXONOMY['RS_GCF_1'] + '\n')
    rdir = os.path.join(directory, 'reads')
    os.mkdir(rdir)
    _write_reads(os.path.join(rdir, 'sample-1'), rng, genomes['GCF_000001.1'], 4000)
    _write_reads(os.path.join(rdir, 'sample-2'), rng,
                 genomes['GCF_000001.1'] + genomes['GCF_000002.1'], 4000)
    _write_reads(os.path.join(rdir, 'sample-3'), rng, other, 4000)
    manifest = os.path.join(directory, 'manifest.tsv')
    with open(manifest, 'w') as fh:
        fh.write('sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n')
        for s in ['sample-1', 'sample-2', 'sample-3']:
            fh.write(f'{s}\t{rdir}/{s}_R1.fastq.gz\t{rdir}/{s}_R2.fastq.gz\n')
    return db, manifest


_SYLPH_DIR = None


def _sylph_data():
    global _SYLPH_DIR
    if _SYLPH_DIR is None:
        _SYLPH_DIR = tempfile.mkdtemp(prefix='la-confluence-sylph-')
        synthetic_sylph_data(_SYLPH_DIR)
    return _SYLPH_DIR


def sylph_database_factory():
    from la_confluence._sylph import SylphDatabaseDirFmt
    directory = _sylph_data()
    db = SylphDatabaseDirFmt()
    for f in ['database.syldb', 'taxonomy.tsv.gz']:
        shutil.copy(os.path.join(directory, 'db', f), os.path.join(str(db), f))
    return Artifact.import_data('SylphDatabase', db)


def paired_reads_factory():
    directory = _sylph_data()
    return Artifact.import_data('SampleData[PairedEndSequencesWithQuality]',
                                os.path.join(directory, 'manifest.tsv'),
                                'PairedEndFastqManifestPhred33V2')


def profile_reads_example(use):
    reads = use.init_artifact('reads', paired_reads_factory)
    database = use.init_artifact('sylph-database', sylph_database_factory)
    use.action(
        use.UsageAction('confluence', 'profile_reads'),
        use.UsageInputs(reads=reads, database=database, threads=1),
        use.UsageOutputNames(taxonomic_abundance='sylph-taxonomic-abundance',
                             sequence_abundance='sylph-sequence-abundance',
                             taxonomy='sylph-taxonomy', stats='sylph-stats'),
    )


# --- paired study: six samples profiled by both sources ----------------------

_PAIRED_SAMPLES = [f'sample-{i}' for i in range(1, 7)]


def _paired_counts(rng, bias):
    """Six samples x three taxa; `bias` tilts the third taxon by source so the
    random forest has a signal."""
    base = rng.integers(50, 500, size=(3, 6)).astype(float)
    base[2] *= bias
    return base


def paired_amplicon_table_factory():
    import biom
    import numpy as np
    rng = np.random.default_rng(1)
    counts = _paired_counts(rng, 0.2)
    counts = np.vstack([counts, rng.integers(1, 20, size=(1, 6))])  # unnamed ASV
    table = biom.Table(counts, ['asv-1', 'asv-2', 'asv-3', 'asv-4'],
                       _PAIRED_SAMPLES)
    return Artifact.import_data('FeatureTable[Frequency]', table)


def paired_amplicon_taxonomy_factory():
    return _taxonomy_factory({
        'asv-1': _TAXONOMY['RS_GCF_1'], 'asv-2': _TAXONOMY['GB_GCA_3'],
        'asv-3': _TAXONOMY['RS_GCF_4'], 'asv-4': _ASV_TAXONOMY['asv-3']})


def paired_metagenome_table_factory():
    import biom
    import numpy as np
    rng = np.random.default_rng(2)
    counts = _paired_counts(rng, 5.0) * 1000
    table = biom.Table(counts, ['GCF_1', 'GCA_3', 'GCF_4'], _PAIRED_SAMPLES)
    return Artifact.import_data('FeatureTable[Frequency]', table)


def paired_metagenome_taxonomy_factory():
    return _taxonomy_factory({'GCF_1': _TAXONOMY['RS_GCF_1'],
                              'GCA_3': _TAXONOMY['GB_GCA_3'],
                              'GCF_4': _TAXONOMY['RS_GCF_4']})


def paired_metadata_factory():
    return _metadata_factory(_PAIRED_SAMPLES,
                             site=['compost'] * 3 + ['soil'] * 3)


def paired_report_example(use):
    tree = use.init_artifact('reference-tree', reference_tree_factory)
    taxonomy = use.init_artifact('reference-taxonomy',
                                 reference_taxonomy_factory)
    genus_tree, = use.action(
        use.UsageAction('confluence', 'gtdb_tree_collapse'),
        use.UsageInputs(tree=tree, taxonomy=taxonomy, level=6),
        use.UsageOutputNames(collapsed_tree='genus-tree'),
    )
    tables = use.construct_artifact_collection('tables', {
        'amplicon': use.init_artifact('amplicon-table',
                                      paired_amplicon_table_factory),
        'metagenome': use.init_artifact('metagenome-table',
                                        paired_metagenome_table_factory)})
    taxonomies = use.construct_artifact_collection('taxonomies', {
        'amplicon': use.init_artifact('amplicon-taxonomy',
                                      paired_amplicon_taxonomy_factory),
        'metagenome': use.init_artifact('metagenome-taxonomy',
                                        paired_metagenome_taxonomy_factory)})
    metadata = use.construct_artifact_collection('metadata', {
        'study': use.init_artifact('study-metadata', paired_metadata_factory)})
    _, frequency_table, _, merged_metadata = use.action(
        use.UsageAction('confluence', 'integrate'),
        use.UsageInputs(tables=tables, taxonomies=taxonomies, tree=genus_tree,
                        level=6, total=1_000_000, metadata=metadata),
        use.UsageOutputNames(relative_table='genus-relative-table',
                             frequency_table='genus-frequency-table',
                             taxonomy='genus-taxonomy',
                             metadata='genus-metadata'),
    )
    use.action(
        use.UsageAction('confluence', 'paired_report'),
        use.UsageInputs(table=frequency_table, tree=genus_tree,
                        metadata=merged_metadata, n_boots=2, permutations=99,
                        cv=3, n_estimators=50, random_state=1),
        use.UsageOutputNames(report='paired-report',
                             distance_matrices='distance-matrices',
                             feature_importance='feature-importance'),
    )
