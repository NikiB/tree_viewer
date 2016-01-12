function Tnew = rotateTree(T,NodeName)
propsForFigure.Name='TmpTree';
fig = view_noa(T,[],propsForFigure)
fig = myview(T,[],propsForFigure);

tr = get(fig,'userdata');

[a b] = ismember(NodeName,tr.names);

hp = b-tr.numLeaves;
ind = 1;
childrenA = false(1,tr.numLabels);
childrenA(tr.tree(hp(ind),1)) = true;
for k = tr.tree(hp(ind),1)-tr.numLeaves:-1:1
    if childrenA(k+tr.numLeaves)
        childrenA(tr.tree(k,:))=true;
    end
end
childrenB = false(1,tr.numLabels);
childrenB(tr.tree(hp(ind),2)) = true;
for k = tr.tree(hp(ind),2)-tr.numLeaves:-1:1
    if childrenB(k+tr.numLeaves)
        childrenB(tr.tree(k,:))=true;
    end
end
permuta = 1:tr.numLabels;
chA = find(childrenA(1:tr.numLeaves));
chB = find(childrenB(1:tr.numLeaves));
if chA(1)<chB(1)
    permuta([chA chB])=[chB chA];
else
    permuta([chB chA])=[chA chB];
end
ipermuta = zeros(1,tr.numLabels);
ipermuta(permuta)=1:tr.numLabels;
tr.names = tr.names(permuta);
tr.dist = tr.dist(permuta);
tr.tree = ipermuta(tr.tree);
tr.par = tr.par(permuta(1:end-1));
tr.selected = tr.selected(permuta);
tr.activeNodes = tr.activeNodes(permuta);
tr.sel2root = tr.sel2root(permuta);
set(fig,'userdata',tr)

activeBranches=find(tr.activeBranches)';
oldPos = tr.y([]); 
    
% propagate last leaf
lastleaf = 1:tr.numLabels;
for ind = tr.numBranches:-1:1
    if ~tr.activeNodes(tr.tree(ind,1))
        lastleaf(tr.tree(ind,:))=lastleaf(ind+tr.numLeaves);
    end
end

% find x coordinates of branches
tr.x = tr.dist; 
for ind = tr.numBranches:-1:1
    tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
end

% find y coordinates of branches
dummy = lastleaf([true,diff(lastleaf(1:tr.numLeaves))~=0]);
tr.y=zeros(tr.numLabels,1);
tr.y(dummy)=1:length(dummy);
for ind = activeBranches
    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
end

% update right labels
todis = tr.names(dummy);
set(tr.ha,'ytick',1:length(dummy),'yticklabel',todis)
tr.maxLabelExtent = max(tr.labelExtents(dummy));

% show only active branches
set(tr.hlines,'Visible','off')
set(tr.hlines(tr.activeNodes),'Visible','on')

% update coordinates in lines
for ind = 1:tr.numLabels-1
     set(tr.hlines(ind),'Ydata',tr.y([ind,ind,tr.par(ind)]))
     set(tr.hlines(ind),'Xdata',tr.x([ind,tr.par([ind ind])]))
end
set(tr.hlines(tr.numLabels),'Ydata',tr.y(tr.numLabels)*[1 1 1])
set(tr.hlines(tr.numLabels),'Xdata',tr.x(tr.numLabels)*[1 1 1])
          
% update dots
mask = false(tr.numLabels,1); mask(1:tr.numLeaves) = true;
set(tr.hdots(1),'Ydata',tr.y(tr.activeNodes&~mask),'Xdata',tr.x(tr.activeNodes&~mask))
set(tr.hdots(2),'Ydata',tr.y(tr.activeNodes&mask),'Xdata',tr.x(tr.activeNodes&mask))            

% update red dots
set(tr.hseldots(1),'Ydata',tr.y(tr.activeNodes&~mask&tr.selected),'Xdata',tr.x(tr.activeNodes&~mask&tr.selected))
set(tr.hseldots(2),'Ydata',tr.y(tr.activeNodes&mask&tr.selected),'Xdata',tr.x(tr.activeNodes&mask&tr.selected))            

% set the axis holders          
set(tr.axhold,'Ydata',[0.5,max(tr.y(tr.activeNodes))+0.5])
if numel(oldPos)
    set(tr.ha,'ylim',get(tr.ha,'ylim')+tr.y([])-oldPos);
else
    set(tr.ha,'ylim',get(tr.ha,'ylim')) % just touch 'YLim' such that the listener is triggered
end

% turn on indicative modes
tr.indicativeMode = false;
set([tr.datatip tr.hpathline],'visible','off')
set(tr.hlines,'color','black')
set(tr.hldots(1),'Xdata',[],'Ydata',[])
set(tr.hldots(2),'Xdata',[],'Ydata',[])

Tnew = phytree(tr.tree,tr.dist,tr.names);
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('TmpTree', names, 15));
close(child_handles(k))