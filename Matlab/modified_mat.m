function [mod_mat, mod_groups] = modified_mat(mat, groups, cell_names, group_dev )


vec = [];
vec_names = [];
temp = 1;

if group_dev
    for i=1:length(groups)
        mod_groups{i} = [temp:temp+length(groups{i})-1];
        temp = temp + length(groups{i});
        %find indices
        cur_group_names = groups{i};
        cur_group_indices = [];
        for j=1:length(cur_group_names)
            cur_name = cur_group_names{j};
            ii = find(strfind(cell_names,cur_name)==1);
            if ~isempty(ii),
                cur_group_indices = [cur_group_indices;ii];
            end            
        end
        
        vec = [vec cur_group_indices];
        groups_title{i} = groups{i}{1};
        groups_title{i} = groups_title{i}(6:end);
    end
    mod_mat = mat(vec,vec);
else
   
end

% Figure
%subplot(2,1,1);
figure;
imagesc(mod_mat); 
colorbar ('WestOutside'); 
axis off;
colors = {'r' [153/255 0 0] 'b' [11/255 132/255 199/255] [0 51/255 102/255] [0 127/255 0] [222/255 125/255 0] 'm' [122/255 16/255 228/255] 'k' 'g' 'y'};
ind_X = 1;
ind_Y = 1;



for i=1:length(groups_inds),
    if group_dev
        cur_names = groups_names{i};
    else
        cur_names = groups_title;
    end
    for j=1:length(cur_names),
        textLabels = text(ind_X,0,cur_names{j},'color',colors{i},'FontSize',8);
        ind_X = ind_X+1;
        set(textLabels(1), 'Rotation',90 );
        textLabels = text(length(vec_names)+1,ind_Y,cur_names{j},'color',colors{i},'FontSize',8);
        ind_Y = ind_Y+1;
    end  
end
if group_dev
    temp = 0;
    for i=1:length(groups_inds)-1
        temp = temp + length(groups_inds{i});
        line([xlim],[temp+0.5 temp+0.5],'LineWidth',2,'Color',[0.7 0.7 0.7]);
        line([temp+0.5 temp+0.5],[ylim],'LineWidth',2,'Color',[0.7 0.7 0.7]);
    end
end
end


% for i=1:length(vec),
%     %text(i,0,num2str(vec(i)));
%     %text(0,i,num2str(vec(i)));
%     text(i,0,(vec_names(i)));
%     text(0,i,(vec_names(i)));
% end
% temp = 0;
% for i=1:length(groups_inds)-1
%     temp = temp + length(groups_inds{i});
%     line([xlim],[temp+0.5 temp+0.5],'LineWidth',4,'Color',[0.7 0.7 0.7]);
%     line([temp+0.5 temp+0.5],[ylim],'LineWidth',4,'Color',[0.7 0.7 0.7]);
% end





% if nargin<5,
%     groups = loci_analysis.group_inds;
% end
% if group_dev   
%     vec = [];
%     temp = 1;
%     for i=1:length(groups)
%         mod_groups{i} = [temp:temp+length(groups{i})-1];
%         temp = temp + length(groups{i});
%         vec = [vec groups{i}];
%     end
%     mat = mat(vec,vec);
% else
%     groups = 1;
% end
% 
% % Figure
% figure;
% imagesc(mat); 
% colorbar ('WestOutside'); 
% %axis off;
% colors = {'r' [153/255 0 0] 'b' [11/255 132/255 199/255] [0 51/255 102/255] [0 127/255 0] [222/255 125/255 0] 'm' [122/255 16/255 228/255] 'k' 'g' 'y'};
% ind_X = 1;
% ind_Y = 1;
% 
% for i=1:length(groups),
%     if group_dev
%         cur_group = groups{i};
%         cur_names = names(cur_group);
%     else
%         cur_group = loci_analysis.group_inds_names;
%         cur_names = names;
%     end
%     for j=1:length(cur_group),
%         textLabels = text(ind_X,0,cur_names{j},'color',colors{i},'FontSize',8);
%         ind_X = ind_X+1;
%         set(textLabels(1), 'Rotation',90 );
%         textLabels = text(length(names)+1,ind_Y,cur_names{j},'color',colors{i},'FontSize',8);
%         ind_Y = ind_Y+1;
%     end  
% end
% if group_dev
%     temp = 0;
%     for i=1:length(groups)-1
%         temp = temp + length(groups{i});
%         line([xlim],[temp+0.5 temp+0.5],'LineWidth',2,'Color',[0.7 0.7 0.7]);
%         line([temp+0.5 temp+0.5],[ylim],'LineWidth',2,'Color',[0.7 0.7 0.7]);
%     end
% end
% end