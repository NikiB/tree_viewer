function [T,groups] = Enrichment_Calculation(T,groups)
%     [QLC_score,lineage_best_scores,ancs_best_score_ind]= Quality_of_Largest_Cluster(T,groups);

    %Entropy_score = Tree_Entropy(T,groups);

    add_internal_pvals = 1;
    [T,node_groups,pvals_groups] = Tree_enrichment_calculation(T,groups,add_internal_pvals);

    groups{end+1} = node_groups;
    groups{end+1} = pvals_groups;

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Save Params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save('Tree.mat',T,'-mat');
    save('Groups.mat',groups,'-mat');