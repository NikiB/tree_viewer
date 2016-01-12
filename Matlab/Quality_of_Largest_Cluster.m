function [QLC_score,lineage_best_scores,ancs_best_score_ind]= Quality_of_Largest_Cluster(T,groups,MIN_LINEAGE_SIZE,MIN_SAME_TYPE_SIZE)

if nargin<4,
    MIN_SAME_TYPE_SIZE = 0.5; %the minimum percent of cells from the same type in the lineage
end
if nargin<3,
    MIN_LINEAGE_SIZE = 0.04; %the minimum percent of cells in the lineage
end


[children,ancs] = get_branch_children(T);
group_amount = length(groups);

%find the ancestors with at least 1% from all leaves
tot_cell_amount = get(T,'NUMLEAVES') - 1;
ancs_indin = [];
for i=1:length(children),
    if length(children{i})>tot_cell_amount*MIN_LINEAGE_SIZE,
        ancs_indin=[ancs_indin i];
    end
end
ancs_amount = length(ancs_indin);

%for each group of cells, find the best lineage
lineage_best_scores = zeros(1,group_amount);
ancs_best_score_ind = zeros(1,group_amount);
for g=1:group_amount,
    cur_group_size = length(groups{g});
    P = zeros(1,ancs_amount);
    Q = zeros(1,ancs_amount);
    for a=1:ancs_amount,
        cur_children = children{ancs_indin(a)};
        inter_group = intersect(cur_children, groups{g});
        P(a) = (length(inter_group) / cur_group_size);
        Q(a) = (length(inter_group) / length(cur_children));
    end
    indin = find(Q >= MIN_SAME_TYPE_SIZE);
    if ~isempty(indin)
        cur_scores = (Q(indin)).*(P(indin));
        [lineage_best_scores(g), ind] = max(cur_scores);
        ancs_best_score_ind(g) = ancs_indin(indin(ind));
    else
        lineage_best_scores(g) = 0;
        ancs_best_score_ind(g) = -1;
    end
end
   
QLC_score = mean(lineage_best_scores);