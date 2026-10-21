# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict

import numpy as np
import pandas as pd
import qiime2
import skbio
from q2_types.tree import NewickFormat

from ._taxonomy import _collapsed_name, _is_named

# GTDB taxonomy strings always have these seven ranks.
GTDB_RANKS = 7


def gtdb_tree_collapse(tree: NewickFormat, taxonomy: pd.Series,
                       level: int) -> (skbio.TreeNode, qiime2.Metadata):
    """The registered action: read the Newick file with underscores kept.

    q2-types' default view converts unquoted underscores to spaces (Newick
    convention), which would turn GTDB tips such as RS_GCF_000006945.2 into
    names the taxonomy does not have. The collapsed names are written quoted,
    so they survive the default reader downstream.
    """
    tree = skbio.TreeNode.read(str(tree), format='newick',
                               convert_underscores=False)
    return _collapse_tree(tree, taxonomy, level)


def _collapse_tree(tree: skbio.TreeNode, taxonomy: pd.Series,
                   level: int) -> (skbio.TreeNode, qiime2.Metadata):
    """Collapse a reference tree to one tip per taxon at `level`.

    Every tip must be in `taxonomy`. Tips are grouped by their taxonomy
    collapsed to `level` (q2-taxa's rule, so the new tip names equal the
    feature IDs of a table collapsed to the same level). A group of one tip
    is renamed; a larger group is sheared off at its most recent common
    ancestor, which must contain no other tips (the group is monophyletic)
    and becomes the tip, keeping its own branch length. The other internal
    node names (GTDB's bootstrap and taxon labels) are kept.

    The log has one row per output tip, keyed by its new name so it serves
    as feature metadata for a table collapsed to the same level: the node's
    original name, the number of genomes it stands for and their IDs.
    """
    if not 1 <= level <= GTDB_RANKS:
        raise ValueError(f'level must be between 1 and {GTDB_RANKS}, '
                         f'got {level}')
    n_ranks = taxonomy.str.count(';') + 1
    if (n_ranks != GTDB_RANKS).any():
        raise ValueError(f'{int((n_ranks != GTDB_RANKS).sum())} taxonomy '
                         f'strings do not have the {GTDB_RANKS} GTDB ranks, '
                         f'e.g. {taxonomy[n_ranks != GTDB_RANKS].iloc[0]!r}')
    tree = tree.copy()
    tips = list(tree.tips())
    missing = [t.name for t in tips if t.name not in taxonomy.index]
    if missing:
        raise ValueError(f'{len(missing)} tree tips are missing from the '
                         f'taxonomy, e.g. {missing[:3]}')
    groups = defaultdict(list)
    for tip in tips:
        groups[_collapsed_name(taxonomy[tip.name], level,
                               GTDB_RANKS)].append(tip)
    unnamed = [n for n in groups if not _is_named(n)]
    if unnamed:
        raise ValueError(f'{len(unnamed)} taxa have no name at level {level}, '
                         f'e.g. {unnamed[:3]}')

    # Find every group's ancestor before editing, on the intact topology.
    plan = []
    for name, nodes in groups.items():
        if len(nodes) == 1:
            plan.append((nodes[0], name, nodes))
            continue
        mrca = tree.lowest_common_ancestor(nodes)
        if mrca.count(tips=True) != len(nodes):
            raise ValueError(f'{name} is not monophyletic: its {len(nodes)} '
                             f'tips share an ancestor with other tips')
        plan.append((mrca, name, nodes))
    log = []
    for node, name, nodes in plan:
        log.append({'feature-id': name,
                    'original-name': node.name if node.name else np.nan,
                    'genomes': len(nodes),
                    'genome-ids': ','.join(n.name for n in nodes)})
        node.children = []
        node.name = name
    # lowest_common_ancestor built skbio's node-lookup cache, which the
    # direct edits above left stale.
    tree.clear_caches()
    log = pd.DataFrame(log).set_index('feature-id').sort_index()
    return tree, qiime2.Metadata(log)
