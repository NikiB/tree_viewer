function [part_set,part_id]=tree_to_partitions_by_name(t,nms)

% B is a numeric array of size [NUMBRANCHES X 2] in which every row represents
% a branch of the tree and it contains two pointers to the branch or leaf nodes, which are its children.
% returnes the binary representation of all partitions obtained by cutting
% internal edges.
% if nms is given, B and leaf_names are reordered. This is for
% bootstrapping in which the randomized resampled trees have different leaf
% ordering - we want a unified one!!!

if nargin<2,
    nms={};
end

B=get(t,'Pointers');
num_BRANCHES=size(B,1);
num_OTU=num_BRANCHES+1;
leaf_names=get(t,'LeafNames');


% reorder according to nms
% if ~isempty(nms),
%     for i=1:length(leaf_names),
%         ind(i)=find(strcmp(nms,leaf_names{i}));
%     end
% end
Bnew=B;
if ~isempty(nms),
    for i=1:size(B,1),
        for j=1:size(B,2),
            if B(i,j)<=num_OTU,
                index=find(strcmp(nms,leaf_names{B(i,j)}));
                Bnew(i,j)=index;
            end
        end
    end
    leaf_names=nms;
    B=Bnew;    
end

num_NODES=num_OTU+num_BRANCHES;

% partitions num_BRANCHES internal branches and for each we have a row vector of 0/1
left_part=zeros(num_NODES,num_OTU);
right_part=zeros(num_NODES,num_OTU);
sons_part=zeros(num_NODES,num_OTU);
left_part(1:num_OTU,1:num_OTU)=eye(num_OTU);
right_part(1:num_OTU,1:num_OTU)=eye(num_OTU);
sons_part(1:num_OTU,1:num_OTU)=eye(num_OTU);
counter=num_OTU+1;
for i=1:num_BRANCHES,
    left_part(counter,:)=sons_part(B(i,1),:);
    right_part(counter,:)=sons_part(B(i,2),:);
    sons_part(counter,:)=left_part(counter,:)+right_part(counter,:);    
    counter=counter+1;
end

% create partition sets
% left_part=left_part(num_OTU+1:end,:);
% right_part=right_part(num_OTU+1:end,:);
% sons_part=sons_part(num_OTU+1:end,:);
% t=2.^(0:num_OTU-1);
% t=repmat(t,num_OTU-1,1);
% part_set=[sum(left_part.*t,2);sum(right_part.*t,2)];
% part_id={};
% counter=1;
% for i=1:size(left_part,1),
%     part_id{counter}=find(left_part(i,:));
%     counter=counter+1;
% end
% for i=1:size(right_part,1),
%     part_id{counter}=find(right_part(i,:));
%     counter=counter+1;
% end

left_part=left_part(num_OTU+1:end,:);
right_part=right_part(num_OTU+1:end,:);
sons_part=sons_part(num_OTU+1:end,:);
% NOA note there is saturation problem in the next line, bootstrap
t=2.^(0:num_OTU-1);
t=repmat(t,num_OTU-1,1);
%part_set=[sum(left_part.*t,2);sum(right_part.*t,2)];
part_set=sum(sons_part.*t,2);
part_id={};
counter=1;
for i=1:size(sons_part,1),
    part_id{counter}=leaf_names(find(sons_part(i,:)));
    counter=counter+1;
end

% for i=1:size(left_part,1),
%     part_id{counter}=find(left_part(i,:));
%     counter=counter+1;
% end
% for i=1:size(right_part,1),
%     part_id{counter}=find(right_part(i,:));
%     counter=counter+1;
% end


