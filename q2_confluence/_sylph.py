# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""sylph profiles as QIIME 2 artifacts.

sylph (Shaw & Yu 2024) sketches reads once per sample and profiles the
sketches by k-mer containment against a database of genomes; sylph-tax
attaches a taxonomy and writes one .sylphmpa file per sample. The parser
below turns those files into a genome-keyed taxonomy and two relative
tables: taxonomic abundance (genome-size corrected, roughly a share of
cells) and sequence abundance (share of reads, renormalized over the
detected genomes when the unknown share was estimated).
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


def _accession(label):
    # TODO: this is AI-generated code; human review pending
    """t__ label (a genome file name) -> accession, e.g. GCF_000227485.1."""
    name = label.split('/')[-1]
    m = re.match(r'^(GC[AF]_\d+\.\d+)', name)
    return m.group(1) if m else re.sub(r'\.(fa|fna|fasta)(\.gz)?$', '', name)


def _read_sylphmpa(path):
    # TODO: this is AI-generated code; human review pending
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
        acc = _accession(strain[3:])
        if acc in rows:
            raise ValueError(f'{path}: accession {acc} appears twice')
        rows[acc] = {'Taxon': ';'.join(ranks), 'rel': float(rel),
                     'seq': float(seq), 'ani': float(ani)}
    return sample, pd.DataFrame.from_dict(rows, orient='index')


def _profiles_to_artifacts(paths, samples=None):
    # TODO: this is AI-generated code; human review pending
    """-> (taxonomic table, sequence table, taxonomy Series, per-sample stats).

    `samples`, when given, is the full sample list: a sample with no
    profile (sylph detected nothing, so sylph-tax wrote no file) gets an
    all-zero column and an unknown share of 1.
    """
    per_sample, taxonomy, stats = {}, {}, []
    for path in sorted(paths):
        sample, df = _read_sylphmpa(path)
        if sample in per_sample:
            raise ValueError(f'sample {sample} appears twice')
        per_sample[sample] = df
        for acc, taxon in df['Taxon'].items():
            if taxonomy.setdefault(acc, taxon) != taxon:
                raise ValueError(f'{acc}: two taxonomy strings across samples')
        stats.append({'sample-id': sample, 'genomes': len(df),
                      'genera': df['Taxon'].str.split(';').str[5].nunique(),
                      'median-ani': df['ani'].median() if len(df) else np.nan,
                      'unknown-share': round(1 - df['seq'].sum() / 100, 6)})
    for sample in (samples or []):
        if sample not in per_sample:
            per_sample[sample] = pd.DataFrame(
                columns=['Taxon', 'rel', 'seq', 'ani'])
            stats.append({'sample-id': sample, 'genomes': 0, 'genera': 0,
                          'median-ani': np.nan, 'unknown-share': 1.0})
    samples, features = sorted(per_sample), sorted(taxonomy)

    def _table(col, renormalize):
        # TODO: this is AI-generated code; human review pending
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
    return _table('rel', False), _table('seq', True), tax, md


# --- database type ----------------------------------------------------------

SylphDatabase = SemanticType('SylphDatabase')


class SylphDbFormat(model.BinaryFileFormat):
    """A sylph genome sketch database (.syldb or the two-stage .syl2db)."""

    def _validate_(self, level):
        # TODO: this is AI-generated code; human review pending
        if os.path.getsize(str(self)) == 0:
            raise model.ValidationError('empty sylph database')


class SylphTaxonomyFormat(model.BinaryFileFormat):
    """sylph-tax's taxonomy metadata: gzipped, two tab-separated columns
    (genome accession, ;-separated taxonomy string)."""

    def _validate_(self, level):
        # TODO: this is AI-generated code; human review pending
        with gzip.open(str(self), 'rt') as fh:
            first = fh.readline()
        if len(first.rstrip('\n').split('\t')) < 2:
            raise model.ValidationError('expected two tab-separated columns')


class SylphDatabaseDirFmt(model.DirectoryFormat):
    database = model.FileCollection(r'database\.syl2?db$', format=SylphDbFormat)
    taxonomy = model.File('taxonomy.tsv.gz', format=SylphTaxonomyFormat)

    @database.set_path_maker
    def _database_path(self, extension):
        # TODO: this is AI-generated code; human review pending
        return f'database.{extension}'

    def database_path(self):
        # TODO: this is AI-generated code; human review pending
        paths = [p for p in os.listdir(str(self))
                 if p in ('database.syldb', 'database.syl2db')]
        if len(paths) != 1:
            raise ValueError(f'expected one database file, found {paths}')
        return os.path.join(str(self), paths[0])


# --- sketches type -----------------------------------------------------------

SylphSketches = SemanticType('SylphSketches')


class SylphSketchFormat(model.BinaryFileFormat):
    """One sample's read sketch (sylph sketch, .sylsp)."""

    def _validate_(self, level):
        # TODO: this is AI-generated code; human review pending
        if os.path.getsize(str(self)) == 0:
            raise model.ValidationError('empty sylph sketch')


class SylphSketchesDirFmt(model.DirectoryFormat):
    """Per-sample read sketches named <sample-id>.paired.sylsp, as sylph
    sketch writes them for paired-end reads named with --lS."""
    sketches = model.FileCollection(r'.+\.paired\.sylsp$',
                                    format=SylphSketchFormat)

    @sketches.set_path_maker
    def _sketch_path(self, sample_id):
        # TODO: this is AI-generated code; human review pending
        return f'{sample_id}.paired.sylsp'

    def sample_paths(self):
        # TODO: this is AI-generated code; human review pending
        """{sample id: sketch path}, sorted by sample id."""
        suffix = '.paired.sylsp'
        return {f[:-len(suffix)]: os.path.join(str(self), f)
                for f in sorted(os.listdir(str(self))) if f.endswith(suffix)}


# --- actions -----------------------------------------------------------------

def _database_url(release):
    # TODO: this is AI-generated code; human review pending
    return ('http://faust.compbio.cs.cmu.edu/sylph-stuff/'
            f'gtdb-r{release}-c200-dbv2.syl2db')


def _taxonomy_url(release):
    # TODO: this is AI-generated code; human review pending
    from sylph_tax.metadata_files import __metadata_file_urls__
    name = f'gtdb_r{release}_metadata.tsv.gz'
    return next(u for u in __metadata_file_urls__ if u.endswith('/' + name))


def _download(url, path):
    # TODO: this is AI-generated code; human review pending
    with urllib.request.urlopen(url) as r, open(path, 'wb') as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)


def fetch_sylph_database(release: str) -> SylphDatabaseDirFmt:
    # TODO: this is AI-generated code; human review pending
    """The pre-built sylph database of the GTDB species representatives
    (one genome per species, sketched at c=200) and sylph-tax's taxonomy
    for the same release."""
    fmt = SylphDatabaseDirFmt()
    _download(_database_url(release), os.path.join(str(fmt), 'database.syl2db'))
    _download(_taxonomy_url(release), os.path.join(str(fmt), 'taxonomy.tsv.gz'))
    return fmt


def sketch_reads(reads: SingleLanePerSamplePairedEndFastqDirFmt,
                 threads: int = 1) -> SylphSketchesDirFmt:
    # TODO: this is AI-generated code; human review pending
    """sylph sketch of every sample's read pairs, named by sample ID."""
    manifest = reads.manifest.view(pd.DataFrame)
    samples = list(manifest.index)
    sketches = SylphSketchesDirFmt()
    with tempfile.TemporaryDirectory() as work:
        lists = {}
        for name, column in [('l1', manifest['forward']),
                             ('l2', manifest['reverse']),
                             ('names', pd.Series(samples, index=samples))]:
            lists[name] = os.path.join(work, name + '.txt')
            with open(lists[name], 'w') as fh:
                fh.write('\n'.join(column.loc[samples]) + '\n')
        _run(['sylph', 'sketch', '-t', str(threads), '--l1', lists['l1'],
              '--l2', lists['l2'], '--lS', lists['names'], '-d', str(sketches)])
    return sketches


def profile_reads(sketches: SylphSketchesDirFmt, database: SylphDatabaseDirFmt,
                  threads: int = 1, estimate_unknown: bool = True
                  ) -> (biom.Table, biom.Table, pd.Series, qiime2.Metadata):
    # TODO: this is AI-generated code; human review pending
    """sylph profile of the sketches, sylph-tax taxprof, then the parser."""
    paths = sketches.sample_paths()
    samples = list(paths)
    with tempfile.TemporaryDirectory() as work:
        profile = os.path.join(work, 'profile.tsv')
        _run(['sylph', 'profile', '-t', str(threads)]
             + (['-u'] if estimate_unknown else [])
             + ['-d', database.database_path(), '-o', profile]
             + list(paths.values()))
        taxprof = os.path.join(work, 'taxprof')
        os.mkdir(taxprof)
        _run(['sylph-tax', 'taxprof', profile, '-t', str(database.taxonomy.view(
              SylphTaxonomyFormat)), '-o', taxprof + '/'],
             env={**os.environ,
                  'SYLPH_TAXONOMY_CONFIG': os.path.join(work, 'config.json')})
        paths = [os.path.join(taxprof, f) for f in os.listdir(taxprof)
                 if f.endswith('.sylphmpa')]
        rel, seq, tax, stats = _profiles_to_artifacts(paths, samples)
    return rel, seq, tax, qiime2.Metadata(stats)


def _read_taxonomy(database):
    # TODO: this is AI-generated code; human review pending
    """{accession: taxonomy string} from the database's sylph-tax taxonomy."""
    path = str(database.taxonomy.view(SylphTaxonomyFormat))
    df = pd.read_csv(path, sep='\t', header=None, index_col=0, dtype=str)
    return df.iloc[:, 0].to_dict()


def _query_to_metadata(path, taxonomy, targets=None):
    # TODO: this is AI-generated code; human review pending
    """sylph query output -> one row per sample and genome, ID
    '<sample>-<accession>' (sample-id is reserved for the ID header), the
    taxonomy joined by accession, and the
    5-95 percentile interval split into two numeric columns (NA -> missing).
    `targets`: keep only rows whose taxonomy or accession contains one."""
    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    if df.empty:
        raise ValueError('sylph query reported no genome with enough sampled '
                         'k-mers in any sample')
    out = pd.DataFrame({'sample': df['Sample_file'],
                        'accession': df['Genome_file'].map(_accession)})
    out['taxonomy'] = out['accession'].map(taxonomy)
    out['adjusted-ani'] = pd.to_numeric(df['Adjusted_ANI'])
    out['effective-coverage'] = pd.to_numeric(df['Eff_cov'])
    interval = df['ANI_5-95_percentile'].str.split('-', n=1, expand=True)
    out['ani-5-percentile'] = pd.to_numeric(interval[0], errors='coerce')
    out['ani-95-percentile'] = pd.to_numeric(interval[1], errors='coerce')
    out['containment-index'] = df['Containment_ind']
    out['naive-ani'] = pd.to_numeric(df['Naive_ANI'])
    if targets:
        lineage = out['taxonomy'].fillna('') + ' ' + out['accession']
        out = out[lineage.apply(lambda s: any(t in s for t in targets))]
        if out.empty:
            raise ValueError(f'no reported genome matches the targets {targets}')
    out = out.sort_values(['sample', 'accession'])
    out.index = pd.Index(out['sample'] + '-' + out['accession'], name='id')
    return qiime2.Metadata(out)


def query_reads(sketches: SylphSketchesDirFmt, database: SylphDatabaseDirFmt,
                threads: int = 1, min_kmers: int = 50,
                targets: list = None) -> qiime2.Metadata:
    # TODO: this is AI-generated code; human review pending
    """sylph query of the sketches against the database, then the parser."""
    with tempfile.TemporaryDirectory() as work:
        out = os.path.join(work, 'query.tsv')
        _run(['sylph', 'query', '-t', str(threads),
              '--min-number-kmers', str(min_kmers),
              '-d', database.database_path(), '-o', out]
             + list(sketches.sample_paths().values()))
        return _query_to_metadata(out, _read_taxonomy(database), targets)


def _run(cmd, env=None):
    # TODO: this is AI-generated code; human review pending
    subprocess.run(cmd, check=True, env=env, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
