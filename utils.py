__author__ = 'veronika'
import json


def load_json(file):
    data = open(file).read()
    json_data = json.loads(data)
    return json_data


def to_hex(n):
    return hex(int(n*255))[2:].upper()


def hex_to_rgb(color_map):
    rgb = '#'
    for i in color_map:
        rgb += to_hex(i)
    return rgb


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


def get_cells_colors(cell_colors_file):
    colors = load_json(cell_colors_file)
    cells_colors = colors['cellColors']
#     cells_colors = strip_keys_to_ids(colors['cellColors'])
    new_dict_colors = dict()
    for n in cells_colors.keys():
        new_dict_colors[n] = hex_to_rgb(cells_colors[n])

    return new_dict_colors


def get_cells_labels(leaf_labels_file):
    labels = load_json(leaf_labels_file)
    cells_labels = labels['leafLabels']
#     cells_labels = strip_keys_to_ids(labels['leafLabels'])
    new_dict_labels = dict()
    for n in cells_labels.keys():
        new_dict_labels[n] = cells_labels[n]

    return new_dict_labels


def get_branch_width(clustering_sizes_file):
    branches = load_json(clustering_sizes_file)
    branches_width = branches['BranchPvals']
    new_dict_width = dict()
    for group in branches_width.keys():
        if not branches_width[group]:
            continue
        for branch in branches_width[group].keys():
            if not branch in new_dict_width.keys():
                new_dict_width[branch] = dict()
            new_dict_width[branch][group] = branches_width[group][branch]

    return new_dict_width


def get_leaf_order(order_vector_file):
    order_vector = load_json(order_vector_file)
    re_order = order_vector['reOrder']
#     re_order = strip_keys_to_ids(order_vector['reOrder'])
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
            if not branch in new_dict_colors.keys():
                new_dict_colors[branch] = dict()
            new_dict_colors[branch][group] = hex_to_rgb(cluster_colors[group][branch])

    return new_dict_colors


def get_legend(legend_file):
    legend_markers = load_json(legend_file)
    legend = legend_markers['Legend']
    new_dict_legend = dict()
    for mark in legend.keys():
        new_dict_legend[mark] = hex_to_rgb(legend[mark])

    return new_dict_legend


def get_bootsrtap_size(bootstrap_sizes_file):
    bootstrap_vector = load_json(bootstrap_sizes_file)
    bootstrap = bootstrap_vector['BranchBootstrap']
    bootstrap_dict = dict()
    for branch in bootstrap.keys():
        if int(bootstrap[branch]) >= 50 and branch != 'Root':
            bootstrap_dict[branch] = int(bootstrap[branch])/10 - 3
        else:
            bootstrap_dict[branch] = 0

    return bootstrap_dict


def get_dup_labels(duplicate_file):
    colors = load_json(duplicate_file)
    duplidate_colors = colors['Duplictaes']
    duplidate_dict_colors = dict()
    for node in duplidate_colors.keys():
        duplidate_dict_colors[node] = hex_to_rgb(duplidate_colors[node])

    return duplidate_dict_colors


def get_description(description_file):
    description_markers = load_json(description_file)
    description = description_markers['Description']
    new_dict_description = dict()
    for mark in description.keys():
        new_dict_description[mark] = description[mark]

    return new_dict_description