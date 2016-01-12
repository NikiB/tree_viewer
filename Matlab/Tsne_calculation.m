function Tsne_calculation(cell_names, mat, Nodes)

%labaling cell names
lbs1 = zeros(length(cell_names),3); % +/-
lbs2 = zeros(length(cell_names),3); % GMP/LMP
lbs3 = zeros(length(cell_names),3); % +/- and GMP/LMP
lbs4 = zeros(length(cell_names),3); % Rel / Diag

[dups_ID,dups_ind,dups_groups] = find_duplicates(cell_names);

for i=1:length(cell_names),
    if ~isempty(strfind(cell_names{i},'-')) | ~isempty(strfind(cell_names{i},'minus'))
        lbs1(i,1:3) = [1 0 0]; % red
        if strfind(cell_names{i},'GMP'),
           lbs3(i,1:3) = [1 0 0]; % red 
        elseif strfind(cell_names{i},'MLP'),
            lbs3(i,1:3) = [0 0 1]; % blue
        end
    elseif ~isempty(strfind(cell_names{i},'+')) | ~isempty(strfind(cell_names{i},'plus')),
        lbs1(i,1:3) = [0 0 1]; % blue
        if strfind(cell_names{i},'GMP'),
            lbs3(i,1:3) = [0 1 0]; % green
        elseif strfind(cell_names{i},'MLP'),
            lbs3(i,1:3) = [0 191/255 191/255]; % turquize
        end
    end
    if strfind(cell_names{i},'GMP'),
        lbs2(i,1:3) = [1 0 0];
    elseif strfind(cell_names{i},'MLP'),
        lbs2(i,1:3) = [0 0 1];
    end
    if strfind(cell_names{i},'Relapse'),
        lbs4(i,1:3) = [1 0 0];
    elseif strfind(cell_names{i},'Diagnosis'),
        lbs4(i,1:3) = [0 0 1];
    end 
end

[mappedX, mapping] = compute_mapping(mat, 'tSNE');


plot_Tsne_map(mappedX,cell_names,Nodes,lbs1,dups_ind,'+ / -');
plot_Tsne_map(mappedX,cell_names,Nodes,lbs2,dups_ind,'GMP / MLP');
plot_Tsne_map(mappedX,cell_names,Nodes,lbs3,dups_ind,'+/- and GMP/MLP');
plot_Tsne_map(mappedX,cell_names,Nodes,lbs4,dups_ind,'Relapse / Diagnosis');

end






function plot_Tsne_map(mappedX,cell_names,Nodes,lbs,dups_ind,tit)

other_marks = {'*' 'square' '>' ' +' 'diamond' '^' 'o'};

figure;
hold on;
draw_tree(mappedX,cell_names,Nodes);
scatter(mappedX(:,1), mappedX(:,2),36,lbs);
for i=1:size(dups_ind,2),
    plot(mappedX(dups_ind(1,i),1), mappedX(dups_ind(1,i),2), 'color',lbs(dups_ind(1,i),:),'marker',other_marks{i},'MarkerSize',10, 'color','k');
    plot(mappedX(dups_ind(2,i),1), mappedX(dups_ind(2,i),2), 'color',lbs(dups_ind(2,i),:),'marker',other_marks{i},'MarkerSize',10, 'color','k' );
    
end
title(tit);

figure;
hold on;
scatter(mappedX(:,1), mappedX(:,2),36,lbs);
for i=1:size(dups_ind,2),
    plot(mappedX(dups_ind(1,i),1), mappedX(dups_ind(1,i),2), 'color',lbs(dups_ind(1,i),:),'marker',other_marks{i},'MarkerSize',10, 'color','k');
    plot(mappedX(dups_ind(2,i),1), mappedX(dups_ind(2,i),2), 'color',lbs(dups_ind(2,i),:),'marker',other_marks{i},'MarkerSize',10, 'color','k' );
    
end
title(tit);

end


