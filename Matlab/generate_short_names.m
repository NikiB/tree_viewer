function [short_cell_names] = generate_short_names(cell_names,cell_data,params)

if ~isempty(params.INFORMATIVE_CELL_NAMES_CRITERION),
    ind_criterion = find(strcmp(cell_data.data_types,params.INFORMATIVE_CELL_NAMES_CRITERION)==1);
    current_vals = cell_data.all_data{ind_criterion};
end

short_cell_names = [];

%ind_cell_ID = find(strcmp(cell_data.data_types,'Cell ID')==1);

for i=1:length(cell_names)
    cur_cell_name = cell_names{i};
    indID = findstr(cur_cell_name,'cID');       
    new_name = [cur_cell_name(indID+3:indID+6)];
    %indID = findstr(cur_cell_name,'Output/');
    %new_name = [(cur_cell_name(indID+7:indID+10))];
    if ~isempty(params.INFORMATIVE_CELL_NAMES_CRITERION)
        ind2 = find(cell_data.cellID == str2num(new_name));
        new_name = [new_name '_' current_vals{ind2}];
    end
    short_cell_names{i} = [new_name '_' cur_cell_name(indID+27:end)];
    
end


[uniques,numUnique] = count_unique(short_cell_names);
ind_dups = find(numUnique>1);
for i=1:length(ind_dups),
    cur_dup_ID = uniques(ind_dups(i));
    ind3 = find(strcmp(short_cell_names,cur_dup_ID)==1);
    for j=1:length(ind3),
        short_cell_names{ind3(j)} = [short_cell_names{ind3(j)} '_dup' num2str(j)];
    end
end



end
