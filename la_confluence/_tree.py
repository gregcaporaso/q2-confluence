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

from ._taxonomy import collapsed_name, is_named


def gtdb_tree_collapse(tree: NewickFormat, taxonomy: pd.Series,
                       level: int) -> skbio.TreeNode:
    """The registered action: read the Newick file with underscores kept.

    q2-types' default view converts unquoted underscores to spaces (Newick
    convention), which would turn GTDB tips such as RS_GCF_000006945.2 into
    names the taxonomy does not have. The collapsed names are written quoted,
    so they survive the default reader downstream.
    """
    tree = skbio.TreeNode.read(str(tree), format='newick',
                               convert_underscores=False, verify=False)
    return collapse_tree(tree, taxonomy, level)


def collapse_tree(tree: skbio.TreeNode, taxonomy: pd.Series,
                  level: int) -> skbio.TreeNode:
    """Collapse a reference tree to one tip per taxon at `level`.

    Every tip must be in `taxonomy`. Tips are grouped by their taxonomy
    collapsed to `level` (q2-taxa's rule, so the new tip names equal the
    feature IDs of a table collapsed to the same level). A group of one tip
    is renamed; a larger group is replaced by its most recent common
    ancestor, which must contain no other tips (the group is monophyletic),
    with a branch length of that ancestor's own length plus the mean
    distance from it to the group's tips. Internal node names are cleared
    and a single-child root is collapsed.
    """
    if level < 1:
        raise ValueError(f'level must be at least 1, got {level}')
    n_ranks = int(taxonomy.str.split(';').str.len().max())
    if level > n_ranks:
        raise ValueError(f'level {level} exceeds the {n_ranks} ranks in the '
                         'taxonomy')
    tree = tree.copy()
    tips = list(tree.tips())
    missing = [t.name for t in tips if t.name not in taxonomy.index]
    if missing:
        raise ValueError(f'{len(missing)} tree tips are missing from the '
                         f'taxonomy, e.g. {missing[:3]}')
    groups = defaultdict(list)
    for tip in tips:
        groups[collapsed_name(taxonomy[tip.name], level, n_ranks)].append(tip)
    unnamed = [n for n in groups if not is_named(n)]
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
        depth = statistics.fmean(n.accumulate_to_ancestor(mrca) for n in nodes)
        plan.append((mrca, name, depth))
    for node, name, depth in plan:
        node.children = []
        node.name = name
        node.length = (node.length or 0.0) + depth
    for node in tree.non_tips(include_self=True):
        node.name = None
    while len(tree.children) == 1:
        child = tree.children[0]
        tree.remove(child)
        child.parent = None
        child.length = None
        tree = child
    tree.clear_caches()
    return tree
