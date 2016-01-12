function [T] = Tree_to_Tsne(T,MS_mat,cell_names)
%% Annotation for tree
    [full_MS_table, nodes_names, Unite_MS_Nodes] = MS_annotation_on_tree(T, MS_mat,cell_names);

    Tsne_calculation(nodes_names, full_MS_table, Unite_MS_Nodes);

    bootstraps = [];
    % %bootstraping
    % [T,bootstraps] = tree_bootstraping(T,D,MS_mat,SNP_mat,cell_names,root_name,Params);
    save('Tree.mat',T,'-mat');
end