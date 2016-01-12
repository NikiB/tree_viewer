function [T,bootstraps] = tree_bootstraping(T,D,MS_mat,SNP_mat,cell_names,root_name,params)

if params.BOOT_ITER==0,
   bootstraps = [];
   return;
end

[part_set_real part_id]=tree_to_partitions_by_name(T,cell_names);


% assign weights to loci 
cW_MS = weight_assignemnt(MS_mat);
cW_SNP = weight_assignemnt(SNP_mat);


% bootstrap values
part_rand = zeros(length(D)-1,params.BOOT_ITER);
ind_all = {};

for i=1:params.BOOT_ITER,
    i
    % resample the data
    [MS_matr] = generate_resample_mat(MS_mat,cW_MS);
    [SNP_matr] = generate_resample_mat(SNP_mat,cW_SNP);
    
    Dr = Distance_calculation(MS_matr,SNP_matr,params);
    
    %reconstruct the tree
    if isempty(params.ROOT),
        [Tr,Z] = my_seqneighjoin(Dr,'equivar',cell_names);
    else
        [Tr,Z] = my_seqneighjoin(Dr,'equivar',cell_names,'reroot',root_name);
    end
         
    %Tr_rand{i} = Tr;
    %Dr_rand{i} = Dr;
    %MS_matr_rand{i} = MS_matr;
    %SNP_matr_rand{i} = SNP_matr;
    [part_rand(:,i),part_id_rand{i}] = tree_to_partitions_by_name(Tr,cell_names);    
end

bootstraps = bootstrap_vals(part_set_real,part_rand);


% % finally change branch names to the bootstrap and Buneman values.
% [children,ancs]=get_branch_children(T);
% 
% 
% % % replace each branch with the bootstrap value
% % leaf_names = get(T,'LeafNames');
% % branch_index = zeros(1,length(children));
% % for i=1:length(part_id),
% %     for j=1:length(children),
% %         if isempty(setdiff(part_id{i},children{j})) & isempty(setdiff(children{j},part_id{i})),
% %             branch_index(j)=i;
% %         end
% %     end
% % end
% 
% node_names = get(T,'NodeNames');
% 
% % branch_names = get(T,'BranchNames');
% % branch_index = 1:length(branch_names);
% % 
% % for branch_i=1:length(branch_names)-1,
% %     ind_branch_i = get(T,'NumLeaves')+branch_i;
% %     br_name = node_names{ind_branch_i};
% %     br_name = br_name(7:end);
% %     br_name = ['br' br_name];
% %     %node_names{ind_branch_i}=[br_name ',bs=' num2str(bootstraps(branch_i))];
% %     node_names{ind_branch_i} = ['bs' br_name '=' num2str(bootstraps(branch_i))];
% % end
% % 
% %T = phytree(get(T,'Pointers'),get(T,'Distances'),node_names);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cW] = weight_assignemnt(mat)

% assign weights to loci so that we don't sample non-informative loci
nonNan_mat = ~isnan(mat);
W = sum(nonNan_mat);
W = W/sum(W);
cW = cumsum(W);

% w=sum(mat>NAN_IDENTIFIER);
% w=w/sum(w);
% cW=cumsum(w);
end


function [matr] = generate_resample_mat(mat,cW)

% resample the data
ind = 1 + floor(size(mat,2)*rand(1,size(mat,2)));

for k = 1:size(mat,2),
    r = rand;
    f = find(r<cW);
    if isempty(f),
        f = 1;
    end
    ind(k) = min(f);
end

matr = mat(:,ind);
end
