__author__ = 'veronika'

import argparse
import ete3
import math
from tree_viewer import utils
import os
import sys
from tree_viewer.cluster_nodes import size_clustering, color_clustering


def node_check(name, t):
    nodes = t.search_nodes(name=name)
    assert len(nodes) >= 1
    if len(nodes) == 0:
        Warning("The json file contains too many nodes")
        return None
    node = nodes[0]
    return node


def tree_draw(tree_file,
              tree_name=None,
              order_vector_file=None,
              cell_colors_file=None,
              clustering_colors_file=None,
              clustering_sizes_file=None,
              intermediate_node_sizes_file=None,
              intermediate_node_labels_file=None,
              leaf_labels_file=None,
              legend_file=None,
              duplicate_file=None,
              tree_scale='linear',
              tree_rotation=True,
              font_size=7,
              font_legend=7,
              node_size=3,
              scale_rate=None,
              distance_factor=1,
              y_scale=True
              ):

    t = ete3.Tree(newick=tree_file, format=1)
    ts = ete3.TreeStyle()
    if tree_rotation:
        ts.rotation = 90
    ts.show_leaf_name = True
    ts.show_scale = False
    ts.scale = 1
    if tree_name:
        ts.title.add_face(ete3.TextFace(tree_name, fsize=20), column=0)

    styles = {}
    max_dist = 0

    # initialize all nodes and branches
    for n in t.traverse():
        styles[n.name] = dict()
        styles[n.name]['style'] = ete3.NodeStyle()
        styles[n.name]['style']['fgcolor'] = 'black'
        styles[n.name]['style']["vt_line_width"] = 1
        styles[n.name]['style']["hz_line_width"] = 1
        max_dist = max(max_dist, n.dist)

    # calculate the scale for the tree (log, linear and right size)
    if tree_scale == 'log':
        max_dist = 0
    for n in t.traverse():
        if tree_scale == 'log':
            root = t.get_tree_root()
            if n == root:
                styles[n.name]['dist'] = 0
            else:
                father_path = 0
                for ancestor in n.get_ancestors():
                    father_path += styles[ancestor.name]['dist']

                dist = math.log10(n.get_distance(root)*distance_factor+1)-father_path
                if dist < 0:
                    dist = 0
                styles[n.name]['dist'] = dist
                max_dist = max(max_dist, dist)

        elif tree_scale == 'linear':
            if max_dist > 1:
                styles[n.name]['dist'] = round(n.dist/max_dist)
            else:
                styles[n.name]['dist'] = n.dist

    # leaf styles and update distance
    if not scale_rate:
        scale_rate = max(1000, round(1/max_dist))
    for n in t.traverse():
        if 'dist' in styles[n.name]:
            n.dist = styles[n.name]['dist']*scale_rate
        if not n.is_leaf():
            styles[n.name]['style']["size"] = 0
        else:
            styles[n.name]['style']["size"] = node_size

    # add bootstrap values to the branches (size of the node)
    if intermediate_node_sizes_file:
        bootsrtap_sizes = utils.get_bootsrtap_size(intermediate_node_sizes_file)
        for branch, size in bootsrtap_sizes.items():
            styles[branch]['style']["size"] = size
            styles[branch]['style']['fgcolor'] = 'black'

    # add colors to the leafs
    if cell_colors_file:
        cells_colors = utils.get_cells_colors(cell_colors_file)
        for name, color in cells_colors.items():
            styles[name]['style']['fgcolor'] = color

    # reorder the tree by pre-proses if possible
    if order_vector_file:
        leaf_order = utils.get_leaf_order(order_vector_file)
        for n in t.traverse('postorder'):
            if n.get_descendants():
                a = ''
                for leaf in n.get_descendants(strategy='postorder'):
                    if leaf.is_leaf():
                        if not a:
                            a = leaf
                b = n.get_descendants(strategy='preorder')[-1]

                if a.is_leaf() and b.is_leaf():
                    if leaf_order[a.name] > leaf_order[b.name]:
                        left, right = n.children
                        n.children = [right, left]

    # add width to branches
    if clustering_sizes_file:
        t, styles = size_clustering(t, styles, clustering_sizes_file)

    # add colors to branches
    if clustering_colors_file:
        t, ts, styles = color_clustering(t, ts, styles, clustering_colors_file)

    # add new leaf labels
    if leaf_labels_file:
        cells_labels = utils.get_cells_labels(leaf_labels_file)
        ts.show_leaf_name = False
        for name, label in cells_labels.items():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            if name in cells_colors:
                name_face = ete3.faces.TextFace(cells_labels[name], fsize=font_size, fgcolor=cells_colors[name])
            else:
                name_face = ete3.faces.TextFace(cells_labels[name], fsize=font_size)

            name_face.margin_left = 3
            node.add_face(name_face, column=0)

    # add duplicate tags to nodes
    if duplicate_file:
        dup_labels = utils.get_dup_labels(duplicate_file)
        for name, color in dup_labels.items():
            node = node_check(name, t)
            if not node:
                continue
            dup_face = ete3.faces.TextFace('*', fsize=10, fgcolor=color)
            dup_face.margin_left = 5
            node.add_face(dup_face, column=1)

    # add y_scale
    if y_scale:
        pass

    # add legend to the tree
    if legend_file:
        legend = utils.get_legend(legend_file)
        for mark in list(legend.keys()):
            ts.legend.add_face(ete3.faces.CircleFace(2, legend[mark]), column=0)
            legend_txt = ete3.faces.TextFace(mark, fsize=font_legend)
            legend_txt.margin_left = 5
            ts.legend.add_face(legend_txt, column=1)
        ts.legend_position = 4

    # set all the styles
    for n in t.traverse():
        if n.name == 'IDroot':
            n.dist = 0
            n.delete()
        if n.is_root():
            n.dist = 0
            n.delete()
        n.set_style(styles[n.name]['style'])
    root = ete3.faces.CircleFace(2, 'white')
    root.border.width = 1
    root.border.color = 'black'
    t.add_face(root, column=0, position='float')

    #t.render("%%inline", tree_style=ts)
    return t, ts


def main():
    parser = argparse.ArgumentParser(description='Draw a tree from JSon files with newick backbone')
    parser.add_argument('-i', '--input', type=str, dest='tree_file', required=True, help='newick tree file distanation')
    parser.add_argument('-o', '--output', type=str, dest='output_file', required=True, help='name for the tree (png file)')
    parser.add_argument('-n', '--name', type=str, dest='tree_name', default=None, help='name for the tree - Title')
    parser.add_argument('-v', '--vectorFile', type=str, dest='order_vector_file', default=None, help='path for order vector file')
    parser.add_argument('-c', '--colorFile', type=str, dest='cell_colors_file', default=None, help='path for cell colors file')
    parser.add_argument('-C', '--clusteringFile', type=str, dest='clustering_colors_file', default=None, help='path for clustering colors file')
    parser.add_argument('-b', '--clusteringSize', type=str, dest='clustering_sizes_file', default=None, help='path for bootstrapping file')
    parser.add_argument('-s', '--nSizeFile', type=str, dest='intermediate_node_sizes_file', default=None, help='path for intermediate node sizes file')
    parser.add_argument('-l', '--nLabelFile', type=str, dest='intermediate_node_labels_file', default=None, help='path for intermediate node labels file')
    parser.add_argument('-L', '--leafFile', type=str, dest='leaf_labels_file', default=None, help='path for leaf labels file')
    parser.add_argument('-D', '--legendFile', type=str, dest='legend_file', default=None, help='path for legend file for the tree')
    parser.add_argument('-d', '--duplicateFile', type=str, dest='duplicate_file', default=None, help='path for duplicates file for the tree')
    parser.add_argument('-S', '--scale', type=str, dest='tree_scale', default='linear', help='choose the scale for the tree linear/log (default=linear)')
    parser.add_argument('-w', '--width', type=int, dest='fig_width', default=None, help='width for the saved figure (default=None)')
    parser.add_argument('-he', '--height', type=int, dest='fig_height', default=None, help='height for the saved figure (default=None)')
    parser.add_argument('-dp', '--dpi', type=int, dest='fig_dpi', default=900, help='DPI resolution for the figure (default=900)')
    parser.add_argument('-r', '--rotation', type=bool, dest='tree_rotation', default=True, help='rotation of the figure (default=True)')
    parser.add_argument('-f', '--fontsize', type=int, dest='font_size', default=7, help='font size for the labels (default=7)')
    parser.add_argument('-fl', '--fontlegend', type=int, dest='font_legend', default=7, help='font size for the legend (default=7)')
    parser.add_argument('-ns', '--nodesize', type=int, dest='node_size', default=3, help='sizes of the leaves (default=3)')
    parser.add_argument('-sr', '--scalerate', type=int, dest='scale_rate', default=None, help='the Y-scale rate (default=None)')
    parser.add_argument('-y', '--yscale', type=str, dest='y_scale', default=True, help='show the Y-scale (default=True)')
    parser.add_argument('-df', '--distancefactor', type=int, dest='distance_factor', default=None, help='the distance factor for small distances(default=1)')


    # tree_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered.newick'
    # cell_colors_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered_LeafColor.json'
    # clustering_colors_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered_ClusterColor.json'
    # clustering_sizes_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered_ClusterWidth.json'
    # leaf_labels_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered_LeafLabel.json'
    # order_vector_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered_LeafOrder.json'
    # intermediate_node_sizes_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_Bootstrap.json'
    # legend_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Ruby/NSR5_AC_X_mat_1a__0_01__30Ruby_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered_Legend.json'
    # output_file = '/net/mraid11/export/data/dcsoft/home/LINEAGE/Hiseq/NSR5/fastq/Calling/Tree_Analysis/Rambam_LCL230/Diag_Rel/oNSR5_AC_X_mat_1a__0_01__30Rambam_transposed_NewNames_filtered_0_0_withRoot_distance_ABS_NJ_reordered_leafTab_fillNAN_1_distance_ABS_NJ_reordered_py.png'
    # tree_name = 'Rubi - fillNaN'

    args = parser.parse_args()
    tree_file = args.tree_file
    output_file = args.output_file
    tree_name = args.tree_name
    order_vector_file = args.order_vector_file
    cell_colors_file = args.cell_colors_file
    clustering_colors_file = args.clustering_colors_file
    clustering_sizes_file = args.clustering_sizes_file
    intermediate_node_sizes_file = args.intermediate_node_sizes_file
    intermediate_node_labels_file = args.intermediate_node_labels_file
    leaf_labels_file = args.leaf_labels_file
    legend_file = args.legend_file
    duplicate_file = args.duplicate_file
    tree_scale = args.tree_scale
    fig_width = args.fig_width
    fig_height = args.fig_height
    fig_dpi = args.fig_dpi
    font_size = args.font_size
    font_legend = args.font_legend
    tree_rotation = args.tree_rotation
    node_size = args.node_size
    scale_rate = args.scale_rate
    y_scale=args.y_scale
    distance_factor = args.distance_factor
    # launch server X
    # sys.setdefaultencoding("utf-8")
    os.environ["DISPLAY"] = ":2"

    tree, ts = tree_draw(tree_file, tree_name=tree_name,
                         order_vector_file=order_vector_file,
                         cell_colors_file=cell_colors_file,
                         clustering_colors_file=clustering_colors_file,
                         clustering_sizes_file=clustering_sizes_file,
                         intermediate_node_sizes_file=intermediate_node_sizes_file,
                         intermediate_node_labels_file=intermediate_node_labels_file,
                         leaf_labels_file=leaf_labels_file,
                         legend_file=legend_file,
                         duplicate_file=duplicate_file,
                         tree_scale=tree_scale,
                         tree_rotation=tree_rotation,
                         font_size=font_size,
                         font_legend=font_legend,
                         node_size=node_size,
                         scale_rate=scale_rate,
                         distance_factor=distance_factor,
                         y_scale=y_scale
                         )

    tree.render(output_file, h=fig_height, w=fig_width, dpi=fig_dpi, tree_style=ts)
    tree.ladderize()
    output_file = output_file.split('.')
    output_file = output_file[0] + '_ladderize.' + output_file[1]
    tree.render(output_file, h=fig_height, w=fig_width, dpi=fig_dpi, tree_style=ts)
    #tree.show(tree_style=ts)
    print('Thank You')

if __name__ == "__main__":
    main()
