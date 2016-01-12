 function [T] = Tree_Plot(T,fig_name,varargin)
    %
    %input: - T: newick tree
    %         fig_name:
    %      varargin - cell_order: vector
    %                 colors: colors by groups
    %                 clustering: cells dict by groups and internal_nodes
    %                 node_size: 
    %                 node_labels: 
    %                 leaf_labels: dictionary by cells ID
    %                 bootstruping: cells dict by groups and internal_nodes
    %                   with bootstrup
    
    if isempty(bootstraps),
        tree_plot_Noa(T,colors,0,'BranchLabels','false','group',groups,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true');
    else
        groups{end+1} = get(T,'BranchNames');
        groups{end+1} = bootstraps;
        tree_plot_Noa_for_bootsraping(T,colors,0,'BranchLabels','false','group',groups,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true');
    end


    %tree_plot_Noa(T,colors,0,'BranchLabels','true','group',groups,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true');
    saveas(gcf,[fig_name '_groups'],'fig');
    generate_tree_legend_group_enrichment(Tree_Params,[fig_name '_groups']);



    %draw tree for duplicate coloring
    [dups_ID,dups_ind,dups_groups] = find_duplicates(get(T,'leafNames'));
    dups_groups{end+1} = [];
    dups_groups{end+1} = [];
    colors = {'b' 'r' 'g' 'm' 'c' 'y'};
    tree_plot_Noa(T,colors,0,'BranchLabels','false','group',dups_groups,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true');
    saveas(gcf,[fig_name '_dups'],'fig');
    generate_tree_legend_PCRdups(Tree_Params,[fig_name '_dups']);
end





function [dups_ID,dups_ind,dups_groups] = find_duplicates(cell_names)
    IDs = [];
    dups_ind = [];
    dups_ID = [];
    dups_groups = [];
    dups_count = 1;

    for i=1:length(cell_names),
        cur_name = cell_names{i};
        inds = strfind(cur_name,'_');
        if ~isempty(inds),
            cur_id = cur_name(1:inds(1)-1);
            ind_dup = find(strcmp(IDs,cur_id)==1);
            if ~isempty(ind_dup),
                dups_ind(1,dups_count) =  ind_dup;
                dups_ind(2,dups_count) =  i;
                dups_ID{dups_count} = cur_id;
                dups_groups{dups_count} = {cell_names{ind_dup} cell_names{i}};
                dups_count = dups_count+1;
            end
            IDs{i} = cur_id;
        end
    end
end