__author__ = 'veronika'

import argparse
import ete2
import math
import utils
import os
import sys


def tree_draw(tree_file,
              tree_name,
              output_file,
              order_vector_file=None,
              cell_colors_file=None,
              clustering_colors_file=None,
              clustering_sizes_file=None,
              intermediate_node_sizes_file=None,
              intermediate_node_labels_file=None,
              leaf_labels_file=None,
              legend_file=None,
              duplicate_file=None
              ):

    t = ete2.Tree(newick=tree_file, format=1)

    ts = ete2.TreeStyle()
    ts.rotation = 90
    ts.show_leaf_name = True
    ts.show_scale = False
    ts.scale = 1
    ts.title.add_face(ete2.TextFace(tree_name, fsize=20), column=0)

    styles = {}
    for n in t.traverse():
        styles[n.name] = dict()
        styles[n.name]['style'] = ete2.NodeStyle()
        styles[n.name]['style']['fgcolor'] = '#FFFFFF'
        if n.dist > 0:
            n.dist = -math.log10(n.dist)*10

        if not n.is_leaf():
            styles[n.name]['style']["size"] = 0

    # add bootstrap values to the branches (size of the node)
    if intermediate_node_sizes_file:
        bootsrtap_sizes = utils.get_bootsrtap_size(intermediate_node_sizes_file)
        for branch, size in bootsrtap_sizes.iteritems():
            styles[branch]['style']["size"] = size
            styles[branch]['style']['fgcolor'] = '#FFFFFF'

    # add colors to the leafs
    if cell_colors_file:
        cells_colors = utils.get_cells_colors(cell_colors_file)
        for name, color in cells_colors.iteritems():
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
                        n.children = [right,left]

    # add width to branches
    if clustering_sizes_file:
        branch_width = utils.get_branch_width(clustering_sizes_file)
        for name, groups in branch_width.iteritems():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            for group, pvalue in groups.iteritems():
                width = min(10,-round(math.log10(pvalue)))
                children = node.get_children()
                for child in children:
                    if group not in styles[child.name]:
                        styles[child.name][group] = dict()
                    styles[child.name][group]["hz_line_width"] = width
                if group not in styles[name]:
                    styles[name][group] = dict()
                styles[name][group]["vt_line_width"] = width
                styles[name]['style']["vt_line_width"] = width

    # add colors to branches
    if clustering_colors_file:
        cluster_colors = utils.get_cluster_colors(clustering_colors_file)
        for name, groups in cluster_colors.iteritems():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            for group, color in groups.iteritems():
                children = node.get_children()
                for child in children:
                    if not group in styles[child.name]:
                        styles[child.name][group] = dict()
                    styles[child.name][group]["hz_line_color"] = color

                    rf = ete2.faces.RectFace(height=styles[child.name][group]["hz_line_width"],
                                             width=child.dist*ts.scale+styles[child.name][group]["hz_line_width"]-0.5,
                                             fgcolor=color,
                                             bgcolor=color)
                    child.add_face(rf,0,position='float')
    #                 styles[child.name][group]["hz_line_width"] = 0
                    styles[child.name]['style']["hz_line_width"] = 0
                styles[name][group]["hz_line_color"] = color
                if not group in styles[name]:
                    styles[name][group] = dict()
                styles[name][group]["vt_line_color"] = color
                styles[name]['style']["vt_line_color"] = color

    # add new leaf labels
    if leaf_labels_file:
        cells_labels = utils.get_cells_labels(leaf_labels_file)
        ts.show_leaf_name = False
        for name, label in cells_labels.iteritems():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            if name in cells_colors:
                name_face = ete2.faces.TextFace(cells_labels[name], fsize=7, fgcolor=cells_colors[name])
            else:
                name_face = ete2.faces.TextFace(cells_labels[name], fsize=7)

            name_face.margin_left = 3
            node.add_face(name_face, column=0)

    if duplicate_file:
        dup_labels = utils.get_dup_labels(duplicate_file)
        for name, color in dup_labels.iteritems():
            nodes = t.search_nodes(name=name)
            assert len(nodes) == 1
            node = nodes[0]
            dup_face = ete2.faces.TextFace('*', fsize=10, fgcolor=color)
            dup_face.margin_left = 5
            node.add_face(dup_face, column=1)

    if legend_file:
        legend = utils.get_legend(legend_file)
        for mark in legend.keys():
            ts.legend.add_face(ete2.faces.CircleFace(2, legend[mark]), column=0)
            legend_txt = ete2.faces.TextFace(mark, fsize=7)
            legend_txt.margin_left = 5
            ts.legend.add_face(legend_txt, column=1)
        ts.legend_position = 4

    for n in t.traverse():
        if n.name == 'L_root':
            n.delete()
        n.set_style(styles[n.name]['style'])
    root = ete2.faces.CircleFace(2, 'white')
    root.border.width = 1
    root.border.color = '#FFFFFF'
    t.add_face(root, column=0, position='float')

    t.ladderize()
    #t.render("%%inline", tree_style=ts)
    return t, ts

def main():
    parser = argparse.ArgumentParser(description='Draw a tree from JSon files with newick backbone')
    parser.add_argument('-i', '--input', type=str, dest='tree_file', help='newick tree file distanation')
    parser.add_argument('-o', '--output', type=str, dest='output_file', help='name for the tree (png file)')
    parser.add_argument('-n', '--name', type=str, dest='tree_name', help='name for the tree - Title')
    parser.add_argument('-v', '--vectorFile', type=str, dest='order_vector_file', default=None, help='path for order vector file')
    parser.add_argument('-c', '--colorFile', type=str, dest='cell_colors_file', default=None, help='path for cell colors file')
    parser.add_argument('-C', '--clusteringFile', type=str, dest='clustering_colors_file', default=None, help='path for clustering colors file')
    parser.add_argument('-b', '--clusteringSize', type=str, dest='clustering_sizes_file', default=None, help='path for bootstrapping file')
    parser.add_argument('-s', '--nSizeFile', type=str, dest='intermediate_node_sizes_file', default=None, help='path for intermediate node sizes file')
    parser.add_argument('-l', '--nLabelFile', type=str, dest='intermediate_node_labels_file', default=None, help='path for intermediate node labels file')
    parser.add_argument('-L', '--leafFile', type=str, dest='leaf_labels_file', default=None, help='path for leaf labels file')
    parser.add_argument('-D', '--legendFile', type=str, dest='legend_file', default=None, help='path for legend file for the tree')
    parser.add_argument('-d', '--duplicateFile', type=str, dest='duplicate_file', default=None, help='path for duplicates file for the tree')


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

    # launch server X
    reload(sys)
    sys.setdefaultencoding("utf-8")
    os.environ["DISPLAY"] = ":2"

    tree, ts = tree_draw(tree_file, tree_name, output_file,
              order_vector_file=order_vector_file,
              cell_colors_file=cell_colors_file,
              clustering_colors_file=clustering_colors_file,
              clustering_sizes_file=clustering_sizes_file,
              intermediate_node_sizes_file=intermediate_node_sizes_file,
              intermediate_node_labels_file=intermediate_node_labels_file,
              leaf_labels_file=leaf_labels_file,
              legend_file=legend_file,
              duplicate_file=duplicate_file
              )

    tree.render(output_file, w=1200, dpi=900, tree_style=ts)
    #tree.show(tree_style=ts)
    print 'Thank You'
