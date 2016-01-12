function [dups_ID,dups_ind,dups_groups] = find_duplicates(cell_names)
IDs = [];
dups_ind = [];
dups_ID = [];
dups_groups = [];
dups_count = 1;

for i=1:length(cell_names),
    
    cur_name = cell_names{i};
    if strfind(cur_name(1:1),'L'),
        inds = strfind(cur_name,'_');
        if ~isempty(inds),
            cur_id = cur_name(inds(1)+1:inds(2)-1);
            ind_dup = find(strcmp(IDs,cur_id)==1);
            if ~isempty(ind_dup),
                dups_ind(1,dups_count) =  ind_dup;
                dups_ind(2,dups_count) =  i;
                dups_ID{dups_count} = cur_id;
                dups_groups{dups_count} = {cell_names{ind_dup} cell_names{i}};
                dups_count = dups_count+1;
            end
            cur_id_num = str2num(cur_id);
            if~isempty(cur_id_num),
                IDs{i} = cur_id;
            end
        end
    end
end

end