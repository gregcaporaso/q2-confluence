# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import statistics
from collections import defaultdict

import pandas as pd
import skbio
from q2_types.tree import NewickFormat

from ._taxonomy import _collapsed_name, _is_named

# GTDB taxonomy strings always have these seven ranks.
GTDB_RANKS = 7


def gtdb_tree_collapse(tree: NewickFormat, taxonomy: pd.Series,
                       level: int) -> skbio.TreeNode:
    # TODO: this is AI-generated code; human review pending
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
                   level: int) -> skbio.TreeNode:
    # TODO: this is AI-generated code; human review pending
    """Collapse a reference tree to one tip per taxon at `level`.

    Every tip must be in `taxonomy`. Tips are grouped by their taxonomy
    collapsed to `level` (q2-taxa's rule, so the new tip names equal the
    feature IDs of a table collapsed to the same level). A group of one tip
    is renamed; a larger group is replaced by its most recent common
    ancestor, which must contain no other tips (the group is monophyletic),
    with a branch length of that ancestor's own length plus the median
    distance from it to the group's tips. Internal node names are cleared.
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

    # Measure on the intact tree, then edit, so distances are not affected
    # by earlier collapses.
    plan = []
    for name, nodes in groups.items():
        if len(nodes) == 1:
            plan.append((nodes[0], name, 0.0))
            continue
        mrca = tree.lowest_common_ancestor(nodes)
        if mrca.count(tips=True) != len(nodes):
            raise ValueError(f'{name} is not monophyletic: its {len(nodes)} '
                             f'tips share an ancestor with other tips')
        depth = statistics.median(n.accumulate_to_ancestor(mrca)
                                  for n in nodes)
        plan.append((mrca, name, depth))
    for node, name, depth in plan:
        node.children = []
        node.name = name
        # The root has no branch length (None); it is the ancestor when a
        # level leaves a single taxon.
        node.length = depth if node.length is None else node.length + depth
    for node in tree.non_tips(include_self=True):
        node.name = None
    tree.clear_caches()
    return tree
