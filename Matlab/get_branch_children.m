function [children,ancs]=get_branch_children(t)

% [children,ancs]=get_branch_children(t)
% 12/2/07
% returnes a cell vector of all ancestor branches for each leaf and all
% leaf children for each branch. Used for bootstrap values

% b = [1 2; 3 4; 5 6; 7 8;9 10];
% t = phytree(b);
B=get(t,'Pointers');
num_leaves=get(t,'NumLeaves');
leaf_names=get(t,'LeafNames');

% for each branch get all the leaf downstream
% first for each leaf keep a vector of all its ancestors up to the root
for i=1:num_leaves,
    ancs{i}=[];
end

for i=1:num_leaves,
    index=i;
    [indexi,indexj]=find(B==index);
    while ~isempty(indexi),
        ancs{i}=[ancs{i} indexi+num_leaves];
        index=indexi+num_leaves;
        [indexi,indexj]=find(B==index);
    end
end

% now find the children of each internal node
for i=num_leaves+1:get(t,'NumNodes'),
    children{i}=[];
end
for i=num_leaves+1:get(t,'NumNodes'),
   for j=1:num_leaves,
       if ~isempty(find(ancs{j}==i)),
           children{i}=[children{i} leaf_names(j)];
       end
   end
end

children=children(num_leaves+1:end);
