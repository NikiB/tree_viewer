function SNP_Analysis(tree, snp_mat,cell_names_SNP,SNP_names,fig_name)

% Input:    tree: phytree
%           snp_mat: raws are cells and cols are anps. Contains numbers
%           which indicates nucleotide: A(1), C(2), G(3), T(4)


leaves_names = get(tree,'LeafNames');
%reorder the snp_mat according to the leaves list
%snp_mat = reorder_SNP_mat(snp_mat,cell_names,leaves_names);

%generate new names to the tree's nodes
nodes_names = get(tree,'NodeNames');
for i=1:length(nodes_names),
    tree_new_names(i) = regexprep(nodes_names(i),' ','_');
end
leaf_amount = get(tree,'numleaves');
for i=1:leaf_amount,
    tree_new_names{i} = ['L_' tree_new_names{i}];
end

tree = phytree(get(tree,'pointers'),get(tree,'distances'),tree_new_names);

%generate new names to the cells with SNPs
for i=1:length(cell_names_SNP),
    new_cell_names_SNP(i) = regexprep(cell_names_SNP(i),' ','_');
    new_cell_names_SNP{i} = ['L_' new_cell_names_SNP{i}];
end
cell_names_SNP = new_cell_names_SNP;

new_cell_names_SNP = intersect(new_cell_names_SNP,tree_new_names);
snp_mat = reorder_SNP_mat(snp_mat,cell_names_SNP,new_cell_names_SNP);

[Nodes,ID] = generate_Node_tree(tree);
[cell_amount, snp_amount] = size(snp_mat);


colors = {'r' 'b' 'g' 'c'};

All_SNP_Trees = [];
All_nodes_names = [];
for i=1:snp_amount,
    
    SNP_Nodes = label_leaves(Nodes,snp_mat(:,i),new_cell_names_SNP);
    
    [cur_SNP_Nodes] = Fitch_Hartigan_algorithm(SNP_Nodes,ID{end});
    All_SNP_Trees{i} = cur_SNP_Nodes;

    %nodes_names = [];
    %[nodes_names] = find_snp_nodes(cur_SNP_Nodes,nodes_names,ID{end},cur_SNP_Nodes.(ID{end}).Label);%     
    %All_nodes_names{i} = nodes_names;
    
    %seperate nodes into group of labels
    
    [groups_for_snp_annotation,groups_leaves] = divied_nodes_to_groups(cur_SNP_Nodes,tree_new_names);
 
    groups_for_snp_annotation{end+1} = [];
    groups_for_snp_annotation{end+1} = [];
    tree_plot_Noa(tree,colors,0,'BranchLabels','false','group',groups_for_snp_annotation,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true');
    title(SNP_names{i});
    cur_fig_name = [fig_name '_SnpAnnot_' SNP_names{i} '_' num2str(i)];
    saveas(gcf,cur_fig_name,'fig');
    saveas(gcf,cur_fig_name,'png');
    
    add_internal_pvals = 1;
    [T,node_groups,pvals_groups] = Tree_enrichment_calculation(tree,groups_leaves,add_internal_pvals);
    groups_leaves{end+1} = node_groups;
    groups_leaves{end+1} = pvals_groups;
    tree_plot_Noa(T,colors,0,'BranchLabels','false','group',groups_leaves,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true');
    title(SNP_names{i});
    cur_fig_name = [fig_name '_SnpEnrichment_' SNP_names{i}  '_' num2str(i)];
    saveas(gcf,cur_fig_name,'fig');
    saveas(gcf,cur_fig_name,'png');
    
    %groups_unite = groups_for_snp_annotation;
    %groups_unite{end-1} = node_groups;
    %groups_unite{end} = pvals_groups;
    %tree_plot_Noa(T,colors,0,'BranchLabels','false','group',groups_unite,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true');

    
    
end

close all;
  
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


function [snp_mat] = reorder_SNP_mat(snp_mat,cell_names_SNP,new_cell_names_SNP)


new_cell_inds = [];
for i=1:length(new_cell_names_SNP),
    ind = find(strcmp(cell_names_SNP,new_cell_names_SNP{i})==1);
    new_cell_inds = [new_cell_inds; ind]
end
snp_mat = snp_mat(new_cell_inds,:);


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



