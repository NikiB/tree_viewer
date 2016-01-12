function draw_tree(mappedX,cell_names,Nodes)

for i=1:length(cell_names),
    cur_location = mappedX(i,:);
    cur_node_name = cell_names{i};
    cur_sons = Nodes.(cur_node_name).Sons;
    if ~isempty(cur_sons),
        for j=1:length(cur_sons),
           son_name =  cur_sons{j};
           son_ind = find(strcmp(cell_names,son_name)==1);
           son_location = mappedX(son_ind,:);
           plot(cur_location(1),cur_location(2),'*k',son_location(1),son_location(2), '*k');
           line([cur_location(1),son_location(1)],[cur_location(2),son_location(2)],'linewidth',0.5,'LineStyle','-','color','k'); 
        end
    end
end

ind_root = find(strcmp(cell_names,'Root')==1);
if ~isempty(ind_root),
    root_location = mappedX(ind_root,:);
    text(root_location(1)+0.01,root_location(2)+0.01 ,'Root');
end

end

