__author__ = 'veronika'


import json
import argparse
import ete2
import math
from Bio import Phylo


def load_tree_with_data(tree_file):
    tree = Phylo.read(tree_file, "newick")
    return tree


def load_json(file):
    data = open(file).read()
    json_data = json.loads(data)
    return json_data


def reorder_tree(order_vector_file):

    reorder_data = load_json(order_vector_file)

    reorder_t = reorder_data
    return reorder_t


def to_hex(n):
    return hex(int(n*255))[2:].upper()


def hex_to_rgb(color_map):
    rgb = '#'
    for i in color_map:
        rgb += to_hex(i)
    return rgb


def get_cells_colors(cell_colors_file):
    colors = load_json(cell_colors_file)
    cells_colors = colors['cellColors']
    new_dict_colors = dict()
    for n in cells_colors.keys():
        new_cell_name = n.split('c')
        new_dict_colors[new_cell_name[1]] = hex_to_rgb(cells_colors[n])

    return new_dict_colors


def get_cluster_colors(clustering_colors_file):
    colors = load_json(clustering_colors_file)
    cluster_colors = colors['BranchColor']
    new_dict_colors = dict()
    for group in cluster_colors.keys():
        if not cluster_colors[group]:
            continue
        for branch in cluster_colors[group].keys():
            new_dict_colors[branch] = hex_to_rgb(cluster_colors[branch])

    return new_dict_colors


def get_branch_width(clustering_sizes_file):
    banches = load_json(clustering_sizes_file)
    banches_width = banches['BranchPvals']
    new_dict_width = dict()
    for group in banches_width.keys():
        if not banches_width[group]:
            continue
        for branch in banches_width[group].keys():
            new_dict_width[branch] = banches_width[group][branch]

    return new_dict_width


def tree_draw(tree_file,
              tree_name,
              order_vector_file,
              cell_colors_file,
              clustering_colors_file,
              clustering_sizes_file,
              intermediate_node_sizes_file,
              intermediate_node_labels_file,
              leaf_labels_file):

    t = ete2.Tree(newick=tree_file, format=1)
    ts = ete2.TreeStyle()
    ts.rotation = 90

    styles = {}
    for n in t.traverse():
        styles[n.name] = ete2.NodeStyle()

    if order_vector_file:
        order = ete2.NodeStyle()
        order = reorder_tree(order_vector_file)

    # add colors to the leafs
    if cell_colors_file:
        cells_colors = get_cells_colors(cell_colors_file)
        for name, color in cells_colors.iteritems():
            styles[name]["fgcolor"] = color

    # add colors to branches
    if clustering_colors_file:
        cluster_colors = get_cluster_colors(clustering_colors_file)
        for name, color in cluster_colors.iteritems():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            children = node.get_children()
            for child in children:
                styles[child.name]["hz_line_color"] = color
            styles[name]["vt_line_color"] = color

    # add width to branches
    if clustering_sizes_file:
        branch_width = get_branch_width(clustering_sizes_file)
        for name, pvalue in branch_width.iteritems():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            width = min(10,-round(math.log10(pvalue)))
            children = node.get_children()
            for child in children:
                styles[child.name]["hz_line_width"] = width
            styles[name]["vt_line_width"] = width

    #
    if intermediate_node_sizes_file:
        pass

    # add
    if intermediate_node_labels_file:
        pass

    # add new leaf labels
    if leaf_labels_file:
        pass

    t.render("%%inline", tree_style=ts)


def analyse_tree(tree_file,
                 tree_name,
                 order_vector_file=None,
                 cell_colors_file=None,
                 clustering_colors_file=None,
                 clustering_sizes_file=None,
                 intermediate_node_sizes_file=None,
                 intermediate_node_labels_file=None,
                 leaf_labels_file=None
                 ):
    """
    preparing all the data for drawing a tree in draw_tree (for now, it's matlab)
    :param tree_file: file name of newick tree
    :param tree_name: name for the tree fig that is generated
    :param order_vector_file:
    :param cell_colors_file:
    :param clustering_colors_file:
    :param intermediate_node_sizes_file:
    :param intermediate_node_labels_file:
    :param leaf_labels_file:
    :param bootstraping_file:
    :return:
    """
    t = load_tree_with_data(tree_file)  # any tree in newick form
    cell_colors = cell_colors_file
    tree_name = ''.format()
    # run matlab groups
    tree_draw(tree_file,
              tree_name,
              order_vector_file,
              cell_colors_file,
              clustering_colors_file,
              clustering_sizes_file,
              intermediate_node_sizes_file,
              intermediate_node_labels_file,
              leaf_labels_file)


if '__main__' == __name__:
    parser = argparse.ArgumentParser(description='Draw a tree from JSon files with newick backbone')
    parser.add_argument('-i', '--input', type=str, dest='tree_file', help='newick tree file distanation')
    parser.add_argument('-n', '--name', type=str, dest='tree_name', help='name for the tree (jpg file)')
    parser.add_argument('-v', '--vectorFile', type=str, dest='order_vector_file', default=None, help='path for order vector file')
    parser.add_argument('-c', '--colorFile', type=str, dest='cell_colors_file', default=None, help='path for cell colors file')
    parser.add_argument('-C', '--clusteringFile', type=str, dest='clustering_colors_file', default=None, help='path for clustering colors file')
    parser.add_argument('-b', '--clusteringSize', type=str, dest='clustering_sizes_file', default=None, help='path for bootstrapping file')
    parser.add_argument('-s', '--nSizeFile', type=str, dest='intermediate_node_sizes_file', default=None, help='path for intermediate node sizes file')
    parser.add_argument('-l', '--nLabelFile', type=str, dest='intermediate_node_labels_file', default=None, help='path for intermediate node labels file')
    parser.add_argument('-L', '--leafFile', type=str, dest='leaf_labels_file', default=None, help='path for leaf labels file')


    args = parser.parse_args()
    tree_file = args.tree_file
    input_file = args.input_file
    tree_name = args.tree_name
    order_vector_file = args.order_vector_file
    cell_colors_file = args.cell_colors_file
    clustering_colors_file = args.clustering_colors_file
    clustering_sizes_file = args.clustering_sizes_file,
    intermediate_node_sizes_file = args.intermediate_node_sizes_file
    intermediate_node_labels_file = args.intermediate_node_labels_file
    leaf_labels_file = args.leaf_labels_file

    analyse_tree(tree_file, tree_name,
                 order_vector_file=order_vector_file,
                 cell_colors_file=cell_colors_file,
                 clustering_colors_file=clustering_colors_file,
                 clustering_sizes_file=clustering_sizes_file,
                 intermediate_node_sizes_file=intermediate_node_sizes_file,
                 intermediate_node_labels_file=intermediate_node_labels_file,
                 leaf_labels_file=leaf_labels_file
                 )

