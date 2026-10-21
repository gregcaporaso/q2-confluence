# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re

import skbio
from qiime2.plugin import model
from q2_types.tree import NewickFormat

# The root of a GTDB tree: a label (quoted or bare) with no branch length,
# right before the closing semicolon, e.g. )'100:d__Bacteria';
_ROOT_LABEL = re.compile(r"\)\s*(?:'[^']*'|[^()\[\]:,;'\s]+)\s*;\s*$")


class GTDBNewickFormat(model.TextFileFormat):
    """A GTDB reference tree as released (e.g. bac120_r232.tree).

    Newick with genome accessions as tips, bootstrap:taxon labels on
    internal nodes, and a quoted root label without a branch length, which
    UniFrac's parser (unifrac-binaries 1.5, run by diversity
    beta-phylogenetic) cannot read. Importing through this format drops that
    one label and changes nothing else; the GTDB file itself is never edited.
    """

    def _validate_(self, level):
        # TODO: this is AI-generated code; human review pending
        with self.open() as fh:
            text = fh.read() if level == 'max' else fh.read(1024)
        if not text.lstrip().startswith('('):
            raise model.ValidationError('not a Newick tree: it does not start '
                                        'with "("')
        if level == 'max':
            if not text.rstrip().endswith(';'):
                raise model.ValidationError('not a Newick tree: it does not '
                                            'end with ";"')
            try:
                skbio.TreeNode.read([text], format='newick',
                                    convert_underscores=False)
            except Exception as e:
                raise model.ValidationError(f'not a Newick tree: {e}')


def _drop_root_label(text: str) -> str:
    # TODO: this is AI-generated code; human review pending
    """The Newick text with a root label that has no branch length removed."""
    return _ROOT_LABEL.sub(');', text.rstrip()) + '\n'


def _gtdb_newick_to_newick(ff: GTDBNewickFormat) -> NewickFormat:
    # TODO: this is AI-generated code; human review pending
    result = NewickFormat()
    with ff.open() as src, result.open() as dst:
        dst.write(_drop_root_label(src.read()))
    return result
