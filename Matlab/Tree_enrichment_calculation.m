function [new_T,node_groups,pvals_groups] = Tree_enrichment_calculation(T,enrichment_groups,add_internal_pvals)

%Calculate p vals
ind_signif = [];
pvalues = [];
for i=1:length(enrichment_groups),
    if ~isempty(enrichment_groups{i})
        % compute hypergeometric enrichment        
            [pvals,alpha_val]=tree_enrichments(T,enrichment_groups{i},0.2);
            if isempty(alpha_val),
                alpha_val=0;
            end
            ind_signif{i} = find(pvals<=alpha_val);
            pvalues{i} = pvals;
        
    end
end




%change branch names to include the hypergeometric pvals
node_names = get(T,'NodeNames');
leaves_num = get(T,'NumLeaves');
% for i=1:length(get(T,'BranchNames')),
%     node_names{leaves_num+i} = '';
% end

for i=1:length(ind_signif),
    cur_pvals = pvalues{i};
    if ~isempty(cur_pvals)
        cur_ind_signif = ind_signif{i};
        if ~isempty(cur_ind_signif)
            for j=1:length(cur_ind_signif),
                ind_branch = get(T,'NumLeaves') + cur_ind_signif(j);
                node_names{ind_branch} = ['G' num2str(i) ',br' num2str(ind_branch) ',P=' num2str(cur_pvals(cur_ind_signif(j)))];
            end
            [most_signif,ind_most_signif] = min(cur_pvals(cur_ind_signif));
            ind_branch = get(T,'NumLeaves') + cur_ind_signif(ind_most_signif);
            node_names{ind_branch} = ['G' num2str(i) ',br' num2str(ind_branch) ',P=' num2str(most_signif)];
        end
    end
end
    
if add_internal_pvals,
    new_T = phytree(get(T,'Pointers'),get(T,'Distances'),node_names);
else
    new_T = T;
end

%For tree painting
for i=1:length(enrichment_groups),
    if isempty(enrichment_groups{i})
        node_groups{i} = [];
        pvals_groups{i} = [];
    else
        node_groups{i} = node_names(get(new_T,'NumLeaves')+ind_signif{i});
        tempvals = pvalues{i};
        pvals_groups{i} = tempvals(ind_signif{i});
    end
end

end
