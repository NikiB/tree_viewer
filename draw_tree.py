__author__ = 'veronika'


import json
import argparse
import ete2
import math


def load_json(file):
    data = open(file).read()
    json_data = json.loads(data)
    return json_data


def reorder_tree(order_vector_file):

    reorder_data = load_json(order_vector_file)

    reorder_t = reorder_data
    return reorder_t


def strip_keys_to_ids(table):
    samp_id = table.keys()[0]
    if samp_id.isdigit():
        return table
    clean_table = dict()
    for key, value in table.iteritems():
        list_key = list(key)
        for i, digit in enumerate(list_key):
            if digit.isdigit():
                new_key = ''.join(list_key[i:])
                break
        clean_table[new_key] = value
    return clean_table


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
        new_dict_colors[n] = hex_to_rgb(cells_colors[n])

    return new_dict_colors


def get_cells_labels(leaf_labels_file):
    labels = load_json(leaf_labels_file)
    cells_labels = strip_keys_to_ids(labels['leafLabels'])
    new_dict_labels = dict()
    for n in cells_labels.keys():
        new_dict_labels[n] = cells_labels[n]

    return new_dict_labels


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


def get_leaf_order(order_vector_file):
    order_vector = load_json(order_vector_file)
    re_order = strip_keys_to_ids(order_vector['reOrder'])
    # re_order = order_vector['reOrder']
    order_dict = dict()
    for n in re_order.keys():
        order_dict[n] = re_order[n]

    return order_dict


def get_cluster_colors(clustering_colors_file):
    colors = load_json(clustering_colors_file)
    cluster_colors = colors['BranchColor']
    new_dict_colors = dict()
    for group in cluster_colors.keys():
        if not cluster_colors[group]:
            continue
        for branch in cluster_colors[group].keys():
            new_dict_colors[branch] = hex_to_rgb(cluster_colors[group][branch])

    return new_dict_colors


def tree_draw(tree_file,
              tree_name,
              order_vector_file=None,
              cell_colors_file=None,
              clustering_colors_file=None,
              clustering_sizes_file=None,
              intermediate_node_sizes_file=None,
              intermediate_node_labels_file=None,
              leaf_labels_file=None
              ):

    t = ete2.Tree(newick=tree_file, format=1)
    ts = ete2.TreeStyle()
    ts.rotation = 90

    styles = {}
    for n in t.traverse():
        styles[n.name] = ete2.NodeStyle()

    # reorder the tree by pre-proses if possible
    if order_vector_file:
        leaf_order = get_leaf_order(order_vector_file)
        for n in t.traverse('postorder'):
            if n.get_descendants():
                a = ''
                for leaf in n.get_descendants(strategy='levelorder'):
                    if leaf.is_leaf():
                        if not a:
                            a = leaf
                b = n.get_descendants(strategy='levelorder')[-1]

                if a.is_leaf() and b.is_leaf():
                    if leaf_order[a.name] > leaf_order[b.name]:
                        left, right = n.children
                        n.children = [right,left]

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
        cells_labels = get_cells_labels(leaf_labels_file)
        ts.show_leaf_name = False
        for name, label in cells_labels.iteritems():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            if cells_colors:
                name_face = ete2.faces.TextFace(cells_labels[name], fsize=8, fgcolor=cells_colors[name])
            else:
                name_face = ete2.faces.TextFace(cells_labels[name], fsize=8)

            node.add_face(name_face, column=0)

    t.render("%%inline", tree_style=ts)


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

    tree_draw(tree_file, tree_name,
                 order_vector_file=order_vector_file,
                 cell_colors_file=cell_colors_file,
                 clustering_colors_file=clustering_colors_file,
                 clustering_sizes_file=clustering_sizes_file,
                 intermediate_node_sizes_file=intermediate_node_sizes_file,
                 intermediate_node_labels_file=intermediate_node_labels_file,
                 leaf_labels_file=leaf_labels_file
                 )

