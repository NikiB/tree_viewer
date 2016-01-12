function [full_MS_table, nodes_names, Unite_MS_Nodes] = MS_annotation_on_tree_ML76(tree, ms_mat,cell_names)

% Input:    tree: phytree
%           


%reorder the snp_mat according to the leaves list

%generate new names to the tree's nodes
nodes_names = get(tree,'NodeNames');
for i=1:length(nodes_names),
    tree_new_names(i) = regexprep(nodes_names(i),' ','_'); 
    tree_new_names(i) = regexprep(tree_new_names(i),'-','_');
end
leaf_amount = get(tree,'numleaves');
for i=1:leaf_amount,
    tree_new_names{i} = ['L_' tree_new_names{i}];
end

tree = phytree(get(tree,'pointers'),get(tree,'distances'),tree_new_names);

%generate new names to the cells with SNPs
for i=1:length(cell_names),
    new_cell_names(i) = regexprep(cell_names(i),' ','_');
    new_cell_names(i) = regexprep(new_cell_names(i),'-','_');

    new_cell_names{i} = ['L_' new_cell_names{i}];
end
cell_names = new_cell_names;

new_cell_names = intersect(new_cell_names,tree_new_names);
ms_mat = reorder_mat(ms_mat,cell_names,new_cell_names);

[Nodes,ID] = generate_Node_tree(tree);
[cell_amount, ms_amount] = size(ms_mat);



Unite_MS_Nodes = Nodes;

for i=1:ms_amount,
    i
    ms_Nodes = label_leaves(Nodes,ms_mat(:,i),new_cell_names);
    
    [cur_ms_Nodes] = Fitch_Hartigan_algorithm(ms_Nodes,ID{end});
    All_SNP_Trees{i} = cur_ms_Nodes;
    
    Unite_MS_Nodes = unite_nodes(Unite_MS_Nodes,cur_ms_Nodes);
    
end

[full_MS_table, nodes_names] = generate_full_MS_table(Unite_MS_Nodes,ms_amount);


  
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nodes,ID] = generate_Node_tree(tree)

[Matrix,ID,dist] = getmatrix(tree);
full_mat = full(Matrix);
Nodes = [];

%start with the root
cell_ind = length(ID);
father_ind = [];
[Nodes] = insert_nodes(Nodes,full_mat,ID,cell_ind,father_ind);

end



function [Nodes] = insert_nodes(Nodes,full_mat,ID,cell_ind,father_ind)
node_name = ID{cell_ind};
if ~isempty(father_ind),
    Nodes.(node_name).Father = ID{father_ind};
else
    Nodes.(node_name).Father = [];
end
Nodes.(node_name).Sons = [];
sons_ind = find(full_mat(cell_ind,:)==1);
if ~isempty(sons_ind),
    Nodes.(node_name).isLeaf = 0;
    current_sons_list = [];
    for i=1:length(sons_ind),
        current_sons_list = [current_sons_list ID(sons_ind(i))];
        [Nodes] = insert_nodes(Nodes,full_mat,ID,sons_ind(i),cell_ind);
    end
    Nodes.(node_name).Sons = current_sons_list;
else
    Nodes.(node_name).isLeaf = 1;
end
Nodes.(node_name).Label = NaN;
end


function [mat] = reorder_mat(mat,cell_names,new_cell_names)


new_cell_inds = [];
for i=1:length(new_cell_names),
    ind = find(strcmp(cell_names,new_cell_names{i})==1);
    new_cell_inds = [new_cell_inds; ind]
end
mat = mat(new_cell_inds,:);


end


function [SNP_Nodes] = label_leaves(Nodes,snp_vec,cell_names_SNP)
SNP_Nodes = Nodes;
for i=1:length(cell_names_SNP),
    SNP_Nodes.(cell_names_SNP{i}).Label = snp_vec(i);
end
end



function [nodes_names] = find_snp_nodes(Nodes,nodes_names,cur_node_name,prev_node_label)

cur_node_label = Nodes.(cur_node_name).Label;

if ~isnan(cur_node_label) & ~isnan(prev_node_label) & cur_node_label~=prev_node_label,
    nodes_names = [nodes_names {cur_node_name}];
end

Sons_list = Nodes.(cur_node_name).Sons;
for i=1:length(Sons_list),
    [nodes_names] = find_snp_nodes(Nodes,nodes_names,Sons_list{i},cur_node_label);
end

end


function [nodes_ind] = find_nodes_inds(ID,nodes_names)

nodes_ind = [];
for i=1:length(nodes_names),
    ind = find(strcmp(ID,nodes_names(i))==1)
    if ~isempty(ind),
       nodes_ind = [nodes_ind ind]; 
    end
end

end


function [groups_all_nodes,groups_leaves] = divied_nodes_to_groups(cur_SNP_Nodes,nodes_names)

groups_all_nodes = [];
groups_leaves = [];
labels_all_nodes = [];
labels_leaves = [];

all_nodes_names = fieldnames(cur_SNP_Nodes);
for i=1:length(all_nodes_names),
   L = cur_SNP_Nodes.(all_nodes_names{i}).Label;
   ind_node = find(strcmp(nodes_names,all_nodes_names{i})==1);
   labels_all_nodes(ind_node) = L;
   if cur_SNP_Nodes.(all_nodes_names{i}).isLeaf,
       labels_leaves(ind_node) = L;
   end
end

for i=1:4,
    indin_all_nodes = find(labels_all_nodes==i);
    groups_all_nodes{i} = nodes_names(indin_all_nodes);
    indin_leaves = find(labels_leaves==i);
    groups_leaves{i} = nodes_names(indin_leaves);    
end

end


function  [Unite_Nodes] = unite_nodes(Unite_Nodes,cur_Nodes),

All_nodes = fieldnames(Unite_Nodes);
for i=1:length(All_nodes),
    new_label = cur_Nodes.(All_nodes{i}).Label;
    Unite_Nodes.(All_nodes{i}).Label = [Unite_Nodes.(All_nodes{i}).Label new_label];
end

end

function [full_MS_table, nodes_names] = generate_full_MS_table(Unite_MS_Nodes,ms_amount)

nodes_names = fieldnames(Unite_MS_Nodes);
full_MS_table = zeros(length(nodes_names),ms_amount);

for i=1:length(nodes_names),
    Labels = Unite_MS_Nodes.(nodes_names{i}).Label;
    full_MS_table(i,:) = Labels(2:end);
end
end



