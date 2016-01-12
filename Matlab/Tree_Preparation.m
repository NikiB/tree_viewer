function [T,MS_mat,SNP_mat,cell_names] = Tree_Preparation(MS_mat,SNP_mat,cell_names,cell_data,Params,fig_name)
% cd('D:/clineage/tree_viewer/test');

% Input:
%       - MS_mat:   MS mutation matrix (cells X loci)
%       - SNP_mat:  SNP mutation matrix (cells X SNPs)
%       - cell_names: list of the cell names, correspoing to MS and SNP mat
%       - Params: tree reconstruction parameter file



    %first,short the names of the cells
    cell_names = generate_short_names(cell_names,cell_data,Params);


    % calculate the root
    [MS_mat,SNP_mat,cell_names,root_name] = root_calculation(MS_mat,SNP_mat,cell_names,Params,cell_data);
    
    %%%%%%%%%%%%%%%%%%%save%%%%%%%%%%%%%%%%%
%     save('root_name.txt',root_name,'-ascii');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %side calculations
    sm = sum(~isnan(MS_mat'));
    fid=fopen([fig_name '_Success_Loci.txt'],'w');
    for i=1:length(sm),
        fprintf(fid,'%s\t%d\n',cell_names{i},sm(i));
    end
    fclose(fid);
    
    %% Find all duplicates
    dups = [];
    counter = 0;
    inds_dup1 = strfind(cell_names,'dup1');
    for i=1:length(inds_dup1),
        if ~isempty(inds_dup1{i}),
            counter = counter+1;
            dup1_cell_name = cell_names{i};
            dup2_cell_name = regexprep(dup1_cell_name,'dup1','dup2');
            ind_dup2 = find(strcmp(cell_names,dup2_cell_name)==1);
            dups(counter,1:2) = [inds_dup1{i},ind_dup2];
        end
    end
    
    %% calculate mu count for dups
    for i=1:size(dups,1),
        diff_sign = MS_mat(dups(i,1),:)-MS_mat(dups(i,1),:);
        mu_count = length(find(~isnan(diff_sign) & diff_sign~=0));
        no_mu_count = length(find(~isnan(diff_sign) & diff_sign==0));
        [mu_count no_mu_count]
    end
    
    %% create mu table for MS table
    mat_mu_count = zeros(size(MS_mat,1),size(MS_mat,1));
    mat_no_mu_count = zeros(size(MS_mat,1),size(MS_mat,1));

    for i=1:size(MS_mat,1),
        for j=i+1:size(MS_mat,1),
            diff_sign = MS_mat(i,:)-MS_mat(j,:);
            mat_mu_count(i,j) = length(find(~isnan(diff_sign) & diff_sign~=0));
            mat_no_mu_count(i,j) = length(find(~isnan(diff_sign) & diff_sign==0));
        end
    end
    mat_mu_count = mat_mu_count + mat_mu_count';
    mat_no_mu_count = mat_no_mu_count + mat_no_mu_count';

    for i=1:size(dups,1),
        [mean(mat_mu_count(dups(i,1),:)) mean(mat_mu_count(dups(i,2),:))]    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %% filter data
    [MS_mat,SNP_mat,cell_names,cell_indin] = filter_data(MS_mat,SNP_mat,cell_names,Params);
    mat_mu_count_after_filtering = mat_mu_count(cell_indin,cell_indin);

    %Tsne_calculation(cell_names, MS_mat)

    %% calculate the distance matrix
    [D] = Distance_calculation(MS_mat,SNP_mat,Params);
    
    %%%%%%%%%%%%%%%%%%%save%%%%%%%%%%%%%%%%%
%     save('DistanceMat.mat',D,'-mat');
%% reconstruct the tree with neighbour joining
    if isempty(Params.ROOT),
        [T,Z] = my_seqneighjoin(D,'equivar',cell_names);
    else
        [T,Z] = my_seqneighjoin(D,'equivar',cell_names,'reroot',root_name);
    end
    
    %save('Tree.mat',T,'-mat');
    phytreewrite('newick_tree_from_mat.newick', T)
    save('Tree_Params.mat',Tree_Params,'-mat');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MS_mat,SNP_mat,leaves_labels,root_name] = root_calculation(MS_mat,SNP_mat,leaves_labels,params,cell_data)

    if isempty(params.ROOT)
        root_name = [];
        return;
    end

    indin = [];
    if ~strcmp(params.ROOT, 'AVE'),
    %     if ~isempty(groups_criterion)
    %         ind_criterion = find(strcmp(cell_data.data_types,groups_criterion)==1);
    %         current_vals = cell_data.all_data{ind_criterion};
    %         for i=1:length(params.ROOT),
    %             indin = [indin ;find(strcmp(current_vals,params.ROOT(i))==1 )];
    %         end
    %     end
    end

    if ~isempty(MS_mat),
        if ~isempty(indin),
            avg_cell_MS = round(nanmean(MS_mat(indin,:),1));
        else
            avg_cell_MS = round(nanmean(MS_mat,1));
        end
        MS_mat = [MS_mat;avg_cell_MS];
    end

    if ~isempty(SNP_mat),
        if ~isempty(indin),
            avg_cell_SNP = round(nanmean(SNP_mat(indin,:),1));
        else
            [vv,avg_cell_SNP] = max(histc(SNP_mat,[1:4]));
        end
        SNP_mat = [SNP_mat;avg_cell_SNP];
    end

    leaves_labels{end+1} = 'Ave_Cell';
    root_name = 'Ave_Cell';
end



function [T,groups] = add_title_to_names(T,groups,groups_title)

    leaves_names = get(T,'leafnames');
    new_leaves_names = leaves_names;
    for i=1:length(groups),
        cur_group = groups{i};
        for j=1:length(cur_group),
            ind = find(strcmp(leaves_names,cur_group{j})==1);
            if ~isempty(ind),
                new_leaves_names{ind} = [new_leaves_names{ind} '_' groups_title{i}];
            end
            cur_group{j} = [cur_group{j} '_' groups_title{i}];
        end
        groups{i} = cur_group;
    end

    T=phytree(get(T,'Pointers'),get(T,'Distances'),new_leaves_names);

end
