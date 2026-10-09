# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""sylph profiles as QIIME 2 artifacts.

sylph (Shaw & Yu 2024) profiles reads by k-mer containment against a
database of genomes; sylph-tax attaches a taxonomy and writes one
.sylphmpa file per sample. The parser below turns those files into a
genome-keyed taxonomy and two relative tables: taxonomic abundance
(genome-size corrected, roughly a share of cells) and sequence abundance
(share of reads, renormalized over the detected genomes when the unknown
share was estimated).
"""

import gzip
import os
import re
import subprocess
import tempfile
import urllib.request

import biom
import numpy as np
import pandas as pd
import qiime2
from qiime2.plugin import SemanticType, model
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt

RANKS = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')


def accession(label):
    """t__ label (a genome file name) -> accession, e.g. GCF_000227485.1."""
    name = label.split('/')[-1]
    m = re.match(r'^(GC[AF]_\d+\.\d+)', name)
    return m.group(1) if m else re.sub(r'\.(fa|fna|fasta)(\.gz)?$', '', name)


def read_sylphmpa(path):
    """-> (sample id, DataFrame indexed by accession: Taxon, rel, seq, ani)."""
    with open(path) as fh:
        header = fh.readline()
    if not header.startswith('#SampleID'):
        raise ValueError(f'{path}: not a .sylphmpa file')
    sample = header.split('\t')[1].strip()
    df = pd.read_csv(path, sep='\t', comment='#', header=0)
    df.columns = [c.split(' ')[0] for c in df.columns]   # 'ANI (if ...)' -> 'ANI'
    rows = {}
    for clade, rel, seq, ani in zip(df['clade_name'], df['relative_abundance'],
                                    df['sequence_abundance'], df['ANI']):
        *ranks, strain = clade.split('|')
        if not strain.startswith('t__'):
            continue
        if len(ranks) != 7 or any(not r.startswith(p) for r, p in zip(ranks, RANKS)):
            raise ValueError(f'{path}: unexpected clade string {clade!r}')
        acc = accession(strain[3:])
        if acc in rows:
            raise ValueError(f'{path}: accession {acc} appears twice')
        rows[acc] = {'Taxon': ';'.join(ranks), 'rel': float(rel),
                     'seq': float(seq), 'ani': float(ani)}
    return sample, pd.DataFrame.from_dict(rows, orient='index')


def profiles_to_artifacts(paths, samples=None):
    """-> (taxonomic table, sequence table, taxonomy Series, per-sample stats).

    `samples`, when given, is the full sample list: a sample with no
    profile (sylph detected nothing, so sylph-tax wrote no file) gets an
    all-zero column and an unknown share of 1.
    """
    per_sample, taxonomy, stats = {}, {}, []
    for path in sorted(paths):
        sample, df = read_sylphmpa(path)
        if sample in per_sample:
            raise ValueError(f'sample {sample} appears twice')
        per_sample[sample] = df
        for acc, taxon in df['Taxon'].items():
            if taxonomy.setdefault(acc, taxon) != taxon:
                raise ValueError(f'{acc}: two taxonomy strings across samples')
        stats.append({'sample-id': sample, 'genomes': len(df),
                      'genera': df['Taxon'].str.split(';').str[5].nunique(),
                      'median_ani': df['ani'].median() if len(df) else np.nan,
                      'unknown_share': round(1 - df['seq'].sum() / 100, 6)})
    for sample in (samples or []):
        if sample not in per_sample:
            per_sample[sample] = pd.DataFrame(
                columns=['Taxon', 'rel', 'seq', 'ani'])
            stats.append({'sample-id': sample, 'genomes': 0, 'genera': 0,
                          'median_ani': np.nan, 'unknown_share': 1.0})
    samples, features = sorted(per_sample), sorted(taxonomy)

    def table(col, renormalize):
        data = np.zeros((len(features), len(samples)))
        for j, s in enumerate(samples):
            v = per_sample[s][col].reindex(features).fillna(0.0).to_numpy() / 100
            if renormalize and v.sum() > 0:
                v = v / v.sum()
            data[:, j] = v
        return biom.Table(data, features, samples)

    tax = pd.Series(taxonomy, name='Taxon').reindex(features)
    tax.index.name = 'Feature ID'
    md = pd.DataFrame(stats).set_index('sample-id').loc[samples]
    return table('rel', False), table('seq', True), tax, md


# --- database type ----------------------------------------------------------

SylphDatabase = SemanticType('SylphDatabase')


class SylphDbFormat(model.BinaryFileFormat):
    """A sylph genome sketch database (.syldb or the two-stage .syl2db)."""

    def _validate_(self, level):
        if os.path.getsize(str(self)) == 0:
            raise model.ValidationError('empty sylph database')


class SylphTaxonomyFormat(model.BinaryFileFormat):
    """sylph-tax's taxonomy metadata: gzipped, two tab-separated columns
    (genome accession, ;-separated taxonomy string)."""

    def _validate_(self, level):
        with gzip.open(str(self), 'rt') as fh:
            first = fh.readline()
        if len(first.rstrip('\n').split('\t')) < 2:
            raise model.ValidationError('expected two tab-separated columns')


class SylphDatabaseDirFmt(model.DirectoryFormat):
    database = model.FileCollection(r'database\.syl2?db$', format=SylphDbFormat)
    taxonomy = model.File('taxonomy.tsv.gz', format=SylphTaxonomyFormat)

    @database.set_path_maker
    def _database_path(self, extension):
        return f'database.{extension}'

    def database_path(self):
        paths = [p for p in os.listdir(str(self))
                 if p in ('database.syldb', 'database.syl2db')]
        if len(paths) != 1:
            raise ValueError(f'expected one database file, found {paths}')
        return os.path.join(str(self), paths[0])


# --- actions -----------------------------------------------------------------

def _database_url(release):
    return ('http://faust.compbio.cs.cmu.edu/sylph-stuff/'
            f'gtdb-r{release}-c200-dbv2.syl2db')


def _taxonomy_url(release):
    from sylph_tax.metadata_files import __metadata_file_urls__
    name = f'gtdb_r{release}_metadata.tsv.gz'
    return next(u for u in __metadata_file_urls__ if u.endswith('/' + name))


def _download(url, path):
    with urllib.request.urlopen(url) as r, open(path, 'wb') as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)


def fetch_sylph_database(release: str) -> SylphDatabaseDirFmt:
    """The pre-built sylph database of the GTDB species representatives
    (one genome per species, sketched at c=200) and sylph-tax's taxonomy
    for the same release."""
    fmt = SylphDatabaseDirFmt()
    _download(_database_url(release), os.path.join(str(fmt), 'database.syl2db'))
    _download(_taxonomy_url(release), os.path.join(str(fmt), 'taxonomy.tsv.gz'))
    return fmt


def profile_reads(reads: SingleLanePerSamplePairedEndFastqDirFmt,
                  database: SylphDatabaseDirFmt, threads: int = 1,
                  estimate_unknown: bool = True
                  ) -> (biom.Table, biom.Table, pd.Series, qiime2.Metadata):
    """sylph sketch, sylph profile, sylph-tax taxprof, then the parser."""
    manifest = reads.manifest.view(pd.DataFrame)
    samples = list(manifest.index)
    with tempfile.TemporaryDirectory() as work:
        lists = {}
        for name, column in [('l1', manifest['forward']),
                             ('l2', manifest['reverse']),
                             ('names', pd.Series(samples, index=samples))]:
            lists[name] = os.path.join(work, name + '.txt')
            with open(lists[name], 'w') as fh:
                fh.write('\n'.join(column.loc[samples]) + '\n')
        sketches = os.path.join(work, 'sketches')
        _run(['sylph', 'sketch', '-t', str(threads), '--l1', lists['l1'],
              '--l2', lists['l2'], '--lS', lists['names'], '-d', sketches])
        profile = os.path.join(work, 'profile.tsv')
        _run(['sylph', 'profile', '-t', str(threads)]
             + (['-u'] if estimate_unknown else [])
             + ['-d', database.database_path(), '-o', profile]
             + sorted(os.path.join(sketches, f) for f in os.listdir(sketches)))
        taxprof = os.path.join(work, 'taxprof')
        os.mkdir(taxprof)
        _run(['sylph-tax', 'taxprof', profile, '-t', str(database.taxonomy.view(
              SylphTaxonomyFormat)), '-o', taxprof + '/'],
             env={**os.environ,
                  'SYLPH_TAXONOMY_CONFIG': os.path.join(work, 'config.json')})
        paths = [os.path.join(taxprof, f) for f in os.listdir(taxprof)
                 if f.endswith('.sylphmpa')]
        rel, seq, tax, stats = profiles_to_artifacts(paths, samples)
    return rel, seq, tax, qiime2.Metadata(stats)


def _run(cmd, env=None):
    subprocess.run(cmd, check=True, env=env, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
