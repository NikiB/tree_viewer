__author__ = 'veronika'

import ete3
import math
from . import utils


def size_clustering(t, styles, clustering_sizes_file):
    # add width to branches
    branch_width = utils.get_branch_width(clustering_sizes_file)
    for name, groups in branch_width.items():
        nodes = t.search_nodes(name=name)
        assert len(nodes) == 1
        node = nodes[0]
        for group, pvalue in groups.items():
            width = min(10, -round(math.log10(pvalue)))
            children = node.get_children()
            for child in children:
                if group not in styles[child.name]:
                    styles[child.name][group] = dict()
                styles[child.name][group]["hz_line_width"] = width
            if group not in styles[name]:
                styles[name][group] = dict()
            styles[name][group]["vt_line_width"] = width
            styles[name]['style']["vt_line_width"] = width
    return t, styles


def color_clustering(t, ts, styles, clustering_colors_file):
    # add colors to branches
    cluster_colors = utils.get_cluster_colors(clustering_colors_file)
    for name, groups in cluster_colors.items():
        nodes = t.search_nodes(name=name)
        assert len(nodes) == 1
        node = nodes[0]
        for group, color in groups.items():
            children = node.get_children()
            for child in children:
                if group not in styles[child.name]:
                    styles[child.name][group] = dict()
                styles[child.name][group]["hz_line_color"] = color

                if 'vt_line_width' in styles[child.name][group]:
                    line_width = child.dist*ts.scale-styles[child.name][group]["vt_line_width"]-0.5
                else:
                    line_width = child.dist*ts.scale
                if 'vt_line_width' in styles[node.name][group]:
                    line_width = line_width + styles[node.name][group]["vt_line_width"]

                rf = ete3.faces.RectFace(height=styles[child.name][group]["hz_line_width"],
                                         width=line_width,
                                         fgcolor=color,
                                         bgcolor=color)
                child.add_face(rf, 0, position='float')
                styles[child.name]['style']["hz_line_width"] = 0
            styles[name][group]["hz_line_color"] = color
            if group not in styles[name]:
                styles[name][group] = dict()
            styles[name][group]["vt_line_color"] = color
            styles[name]['style']["vt_line_color"] = color
    return t, ts, styles
