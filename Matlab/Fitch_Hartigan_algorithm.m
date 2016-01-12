function [ Nodes,label ] = Fitch_Hartigan_algorithm(Nodes,root_name)

% Small Parsimony: Fitch’s Algorithm
% Traverse tree “up”, from leaves to root, finding sets of possible ancestral states (labels) for each internal node.
% Traverse tree “down”, from root to leaves, determining ancestral states (labels) for internal nodes.
% Key observation: Different sites are independent. Can solve one site at a time.

% Bottom-Up phase
Nodes_phase1 = @compute_internal_S1S2_labels;
[Nodes_phase1,label,Nodes_phase1.(root_name)] = compute_internal_S1S2_labels(Nodes,Nodes.(root_name));
Nodes = Nodes_phase1;

% Top-Down phase
Nodes_phase2 = @compute_internal_S1S2_labels;
ancestral_label = [];
[Nodes_phase2,Nodes_phase2.(root_name)] = compute_final_internal_labels(Nodes,Nodes.(root_name),ancestral_label);
Nodes = Nodes_phase2;




end



function [Nodes,label_s1,Node] = compute_internal_S1S2_labels(Nodes,Node)

if Node.isLeaf,
    label_s1 = Node.Label;
else
   labels = [];
   for i=1:length(Node.Sons),
       if any(strcmp(Node.Sons{i},fieldnames(Nodes))),
        [Nodes,label_s1,Nodes.(Node.Sons{i})] = compute_internal_S1S2_labels(Nodes,Nodes.(Node.Sons{i}));
        labels = [labels label_s1];
      end
   end
   %determine the label of the node
   labels = labels(~isnan(labels));
   if isempty(labels),
       labels(1) = NaN;
   end
   all_possibles_labels = unique(labels);
   counts = zeros(1,length(all_possibles_labels));
   for i=1:length(counts),
       counts(i) = length(find(labels==all_possibles_labels(i)));
   end
   max_count = max(counts);
   S1 = all_possibles_labels(find(counts==max_count));
   S2 = all_possibles_labels(find(counts==max_count-1));
   Node.Label_S1 = S1;
   Node.Label_S2 = S2;
   label_s1 = Node.Label_S1;   
end

end


function [Nodes,Node] = compute_final_internal_labels(Nodes,Node,ancestral_label)

if ~Node.isLeaf,
    possible_labels = [];
    if isempty(ancestral_label),
        possible_labels = [Node.Label_S1];
    else
        if ~isempty(find(Node.Label_S1==ancestral_label)),
            possible_labels = ancestral_label;
        elseif ~isempty(find(Node.Label_S2==ancestral_label)),
            possible_labels = [Node.Label_S1,ancestral_label];
        else
            possible_labels = [Node.Label_S1];
        end
    end
%     if length(find(isnan(possible_labels)==1))==length(possible_labels), %%%%
%         Node.Label = ancestral_label;
%     else
        ind = randi(length(possible_labels),1);
        Node.Label = possible_labels(ind);
    %end
    for i=1:length(Node.Sons),
        if any(strcmp(Node.Sons{i},fieldnames(Nodes))),
            [Nodes,Nodes.(Node.Sons{i})] = compute_final_internal_labels(Nodes,Nodes.(Node.Sons{i}),Node.Label);
        end
    end
end
end

