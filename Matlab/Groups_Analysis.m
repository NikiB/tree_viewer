function [groups] = Groups_Analysis(T,cell_data_file_name)
    %% create groups and clusters
    cell_data = read_cell_data_file(cell_data_file_name);
    [groups,groups_cellInds,colors] = groups_division(T,cell_data);
    
    for i=1:length(colors),
        colors{i} = str2num(colors{i});
    end
    save('Groups.mat',groups,'-mat');
    save('Colors.mat',colors,'-mat');
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [groups_cellNames,groups_cellInds,uniques_colors] = groups_division(T,cell_data)
    
    ind_color = find(strcmp(cell_data.data_types, 'Group Color')==1);
    color_vals = cell_data.all_data{ind_color};

    [uniques_colors] = count_unique(color_vals);
    for i=1:length(uniques_colors),
        groups_cellNames{i} = [];
        groups_cellInds{i} = [];
    end

    leaf_names = get(T,'LeafNames');

    for i=1:length(leaf_names),
        cur_name = leaf_names{i};
        ind1 = strfind(cur_name,'_');
        if ~isempty(ind1)
            cur_name = [cur_name(1:ind1(1)-1)];
        end
        n = str2num(cur_name);
        if ~isempty(n)
            ind2 = find(cell_data.cellID==n);
            if ~isempty(ind2)
                cur_color = color_vals(ind2);
                color_group_ind = find(strcmp(uniques_colors,cur_color)==1);
                groups_cellNames{color_group_ind} = [groups_cellNames{color_group_ind} leaf_names(i)];
                groups_cellInds{color_group_ind} = [groups_cellInds{color_group_ind} i];
            end
        end
    end
end
%%

function [data] = read_cell_data_file(cell_data_file_name)

    [num char raw] = xlsread(cell_data_file_name);
    data.cellID = num;

    data.data_types = raw(1,:);

    for i=1:length(data.data_types),
         data.all_data{i} = [];
    end

    for i=1:length(data.data_types),
         data.all_data{i} = raw(2:end,i);
    end

end



