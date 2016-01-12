__author__ = 'veronika'

import scipy.io as sio
import os
from ete2 import Tree, TreeStyle
from subprocess import call
from pymatbridge import Matlab
from utils import user_cells_report


t = Tree("((a,b),c);")

ts = TreeStyle()

# rotate tree
ts.show_leaf_name = True
ts.rotation = 90

# show node names and branch information
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.show(tree_style=ts)



t.render("mytree.png")


def create_json_from_mat(mat, json_name):
    fields = sio.whosmat(mat)
    mat_contents = sio.loadmat(mat)
    data_name = fields[0][0]
    oct_a = mat_contents[data_name]
    with open(json_name, 'w') as f:
        json.dump(oct_a, f)


def create_mat_from_json(json, mat_name):

    json_data = json.loads(json)
    sio.savemat(mat_name, json_data)


def matlab_run():
    matlabFunc = '/clineage/tree_viewer/Matlab/drawTree.m'
    user_cells_report.print_cells_table(partner_name)
    json_file = 'cells_color.JSon'
    create_mat_from_json(json_file, mat_name)

    mlab = Matlab(matlab='/Applications/MATLAB_R2011a.app/bin/matlab')
    mlab.start()
    mlab.run(matlabFunc, {'T': tree_file,
                          'tree_name': tree_name,
                          'order_vector': order_vector_file,
                          'cell_colors_file': cell_colors_file,
                          'clustering_colors_file': clustering_colors_file,
                          'clustering_sizes_file': clustering_sizes_file,
                          'intermediate_node_sizes_file': intermediate_node_sizes_file,
                          'intermediate_node_labels_file': intermediate_node_labels_file,
                          'leaf_labels_file': leaf_labels_file})

    # DosCmd = 'matlab -wait -automation -nosplash -r "run \'' + to_run + "', exit\""
    # # call(['matlab', '-wait', '-automation', '-nosplash', '-r', 'run \' + to_run + \', exit'])
    # os.system(DosCmd)
