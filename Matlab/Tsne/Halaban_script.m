
close all;

addpath('S:\Noa\Tree_analysis_Feb2015\Ver_2.3\');
addpath('S:\Noa\drtoolbox\');


load S:\Noa\Tree_analysis_Feb2015\Ver_2.3\Tsne\Halaban_YUCLAT_tree_data

[full_MS_table, nodes_names, Unite_MS_Nodes] = MS_annotation_on_tree(T, MS_mat, cell_names);

mat = MS_mat;




lbs_cell_names = zeros(length(cell_names),3);
for i=1:length(cell_names),
    cur_cell_name = cell_names{i};
    for j=1:length(groups),
        ind = find(strcmp(groups{j},cur_cell_name)==1);
        if ~isempty(ind),
            lbs_cell_names(i,1:3) = str2num(cell2mat(colors(j)));
        end
    end   
end


for i=1:length(nodes_names),
   temp_nodes_names{i} = regexprep(nodes_names{i},'L_','');
end

lbs_nodes_names = zeros(length(nodes_names),3);
for i=1:length(temp_nodes_names),
    cur_node_name = temp_nodes_names{i};
    for j=1:length(groups),
        ind = find(strcmp(groups{j},cur_node_name)==1);
        if ~isempty(ind),
            lbs_nodes_names(i,1:3) = str2num(cell2mat(colors(j)));
        end
    end
end




for i=1:3,
figure;
[mappedX, mapping] = compute_mapping(mat, 'tSNE',2,length(cell_names));
scatter(mappedX(:,1), mappedX(:,2),36,lbs_cell_names);

figure;
[mappedX, mapping] = compute_mapping(full_MS_table, 'tSNE',2,length(nodes_names));
hold on;
draw_tree(mappedX,nodes_names,Unite_MS_Nodes);
scatter(mappedX(:,1), mappedX(:,2),36,lbs_nodes_names);
end

figure;
[mappedX, mapping] = compute_mapping(mat, 'PCA');
scatter(mappedX(:,1), mappedX(:,2),36,lbs_cell_names);

figure;
[mappedX, mapping] = compute_mapping(full_MS_table, 'PCA');
hold on;
draw_tree(mappedX,nodes_names,Unite_MS_Nodes);
scatter(mappedX(:,1), mappedX(:,2),36,lbs_nodes_names);
