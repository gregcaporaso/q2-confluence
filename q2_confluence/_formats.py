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
# A bare (unquoted) label containing an underscore: a tip after '(' or ','
# or an internal node after ')', up to its branch length or the next
# delimiter; e.g. RS_GCF_000006945.2.
_BARE_UNDERSCORE_LABEL = re.compile(
    r"(?<=[(,)])([^(),:;'\s]*_[^(),:;'\s]*)(?=[:,);])")


class GTDBNewickFormat(model.TextFileFormat):
    """A GTDB reference tree as released: the bacterial bac120_r<release>.tree
    or the archaeal ar53_r<release>.tree (both checked on r226 and r232).

    Newick with genome accessions as tips, bootstrap:taxon labels on
    internal nodes, and a quoted root label without a branch length
    ('100:d__Bacteria', '100:d__Archaea'), which UniFrac's parser
    (unifrac-binaries 1.5, run by diversity beta-phylogenetic) cannot read.
    The artifact imported through this format holds the tree without that
    one label and with every bare label that contains an underscore quoted
    (the tips, RS_GCF_... and GB_GCA_...), so that readers following the
    Newick convention, q2-types' skbio.TreeNode view among them, keep the
    underscores instead of turning them into spaces; otherwise the tree is
    as GTDB wrote it.

    Validation always parses the whole file (2 s for the 9 MB bacterial
    tree), whatever the level.
    """

    def _validate_(self, level):
        with self.open() as fh:
            text = fh.read()
        if not text.lstrip().startswith('('):
            raise model.ValidationError('not a Newick tree: it does not start '
                                        'with "("')
        if not text.rstrip().endswith(';'):
            raise model.ValidationError('not a Newick tree: it does not end '
                                        'with ";"')
        try:
            skbio.TreeNode.read([text], format='newick',
                                convert_underscores=False)
        except Exception as e:
            raise model.ValidationError(f'not a Newick tree: {e}')


def _drop_root_label(text: str) -> str:
    """The Newick text with a root label that has no branch length removed."""
    return _ROOT_LABEL.sub(');', text.rstrip()) + '\n'


def _quote_underscore_labels(text: str) -> str:
    # TODO: this is AI-generated code; human review pending
    """The Newick text with every bare label that contains an underscore
    quoted; labels already quoted, and labels without an underscore, are
    left as they are."""
    return _BARE_UNDERSCORE_LABEL.sub(r"'\1'", text)


def _gtdb_newick_to_newick(ff: GTDBNewickFormat) -> NewickFormat:
    # TODO: this is AI-generated code; human review pending (quoting added on
    # 2026-10-23, after the function had passed)
    result = NewickFormat()
    with ff.open() as src, result.open() as dst:
        dst.write(_quote_underscore_labels(_drop_root_label(src.read())))
    return result
