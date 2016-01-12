close all;

addpath('S:\Noa\Tree_analysis_Feb2015\Ver_2.3\');
addpath('S:\Noa\drtoolbox\');


load S:\Noa\Tree_analysis_Feb2015\Ver_2.3\Tsne\ML61_tree_data

mat_real(find(mat_real<0)) = NaN;

[full_MS_table, nodes_names, Unite_MS_Nodes] = MS_annotation_on_tree(tr_real, mat_real,labels_real);

cell_names = labels_real;
for i=1:length(cell_names),
    new_cell_names(i) = regexprep(cell_names(i),' ','_');
    new_cell_names(i) = regexprep(new_cell_names(i),'-','_');
    new_cell_names{i} = ['L_' new_cell_names{i}];
end
cell_names = new_cell_names;

mat = mat_real;


lbs_cell_names = zeros(length(cell_names),3);
for i=1:length(cell_names),
    if strfind(cell_names{i},'CR'),
        lbs_cell_names(i,1:3) = [0 0 1]; % blue   
    elseif strfind(cell_names{i},'LI3_5'),
        lbs_cell_names(i,1:3) = [1 0 0]; % red
    elseif strfind(cell_names{i},'LI3_6'),
        lbs_cell_names(i,1:3) = [0 204/255 0]; %green
    end
end

lbs_nodes_names = zeros(length(nodes_names),3);
for i=1:length(nodes_names),
    if strfind(nodes_names{i},'CR'),
        lbs_nodes_names(i,1:3) = [0 0 1]; % blue   
    elseif strfind(nodes_names{i},'LI3_5'),
        lbs_nodes_names(i,1:3) = [1 0 0]; % red
    elseif strfind(nodes_names{i},'LI3_6'),
        lbs_nodes_names(i,1:3) = [0 204/255 0]; %green
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
