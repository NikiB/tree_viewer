function [handles,x,y] = python_tree_plot(tree,colors,scale,varargin)
% tree_plot_Noa(T,colors,0,);
% vargins = 'BranchLabels','false','group',groups,'display_branch_labels',1,'ORIENTATION','bottom','TERMINALLABELS','true'

    LINSCALE = 0;
    LOGSCALE = 1;
    MARKER_SIZE=5;
    
%%
%PLOT renders a phylogenetic tree.
%
%   PLOT(TREE) renders a phylogenetic tree object into a MATLAB figure as a
%   phylogram. The significant distances between branches and nodes are in
%   horizontal direction, vertical coordinates are accommodated only for
%   display purposes. Handles to graph elements are stored in the
%   'UserData' figure field, such that graphic properties can be easily
%   modified.
%
%   PLOT(TREE,ACTIVEBRANCHES) hides the non'active branches and all their
%   descendants. ACTIVEBRANCHES is a logical array of size
%   [numBranches x 1] indicating the active branches.
%
%   PLOT(...,'TYPE',type) selects the method to render the phylogenetic
%   tree. Options are: 'square' (default), 'angular', and 'radial'.
%
%   PLOT(...,'ORIENTATION',orient) will orient the phylogenetic tree within
%   the figure window. Options are: 'top', 'bottom', 'left' (default), and,
%   'right'. Orientation parameter is valid only for phylograms or
%   cladograms.
%
%   PLOT(...,'BRANCHLABELS',value) hides/unhides branch labels. Options are
%   true or false. Branch labels are placed next to the branch node.
%   Defaults to false (true) when TYPE is (is not) 'radial'.
%
%   PLOT(...,'LEAFLABELS',value) hides/unhides leaf labels. Options are
%   true or false. Leaf labels are placed next to the leaf nodes. Defaults
%   to false (true) when TYPE is (is not) 'radial'.
%
%   PLOT(...,'TERMINALLABELS',value) hides/unhides terminal labels. Options
%   are true (default) or false. Terminal labels are placed over the axis
%   tick labels, ignored when 'radial' type is used.
%
%   H = PLOT(...) returns a structure with handles to the graph elements.
%
%   Example:
%
%       tr = phytreeread('pf00002.tree')
%       plot(tr,'type','radial')
%
%       % Graph element properties can be modified as follows:
%
%       h=get(gcf,'UserData')
%       set(h.branchNodeLabels,'FontSize',6,'Color',[.5 .5 .5])
%
%   See also PHYTREE, PHYTREE/VIEW, PHYTREEREAD, PHYTREETOOL, SEQLINKAGE.

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.6.10 $ $Author: batserve $ $Date: 2006/06/16 20:06:45 $

    %%
    terminal_colors = colors;
    internal_color = colors;
    if numel(tree)~=1
        error('Bioinfo:phytree:plot:NoMultielementArrays',...
            'Phylogenetic tree must be an 1-by-1 object.');
    end


    % set defaults
    display_branch_labels = 'off';
    display_leaf_labels = true;
    LeafOrder = NaN;
    cellColorsVector = NaN;
    clusterColorsVector = NaN;
    clusterSizeVector = NaN;
    intermediateNodeSize = NaN;
    intermediatoNodeLabel = NaN;
    leafLabels = NaN;
    branchLabels = NaN;
    dispBranchLabels = NaN;
    dispLeafLabels = NaN;
    
    renderType = 'square';
    orientation = 'left';
    rotation = 0;

    %%
    tree = struct(tree);
    tree.numBranches = size(tree.tree,1);

    tree_names2=tree.names;
    % modify branch names for easier visualization
    for i=1:length(tree.names),
        if ~isempty(strfind(tree.names{i},'en=')),
            index=strfind(tree.names{i},'=');
            tree.names{i}=tree.names{i}(index+1:end);
        end
    end

    if nargin>1 && islogical(varargin{1})
        activeBranches = varargin{1};
        argStart = 2;
    else
        activeBranches = true(tree.numBranches,1);
        argStart = 1; 
    end

    %%
    %check all the varargin
    if nargin - argStart > 0
        if rem(nargin - argStart-1,2) == 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%301111%%%%%%%%%%%%%%
            error('Bioinfo:phytree:plot:IncorrectNumberOfArguments',...
                'Incorrect number of arguments to %s.',mfilename);
        end
        okargs = {'type','orientation','rotation',...
            'LeafOrder','cellColor','clusterColor','clusterSize',...
            'intermediateNodeSize','intermediatoNodeLabel','LeafLabels',...
            'branchlabels','display_leaf_labels','Group','display_branch_labels'};
        for j = argStart:2:nargin-argStart-2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%301111%%%%%%%%%%%%%%
            pname = varargin{j};
            pval = varargin{j+1};
            k = find(strncmpi(pname,okargs,numel(pname)));
            if isempty(k)
                error('Bioinfo:phytree:plot:UnknownParameterName',...
                    'Unknown parameter name: %s.',pname);
            elseif length(k)>1
                error('Bioinfo:phytree:plot:AmbiguousParameterName',...
                    'Ambiguous parameter name: %s.',pname);
            else
                switch(k)
                    case 1 % type
                        oktypes={'square','angular','radial'};
                        l = strmatch(lower(pval),oktypes); %#ok
                        if isempty(l)
                            error('Bioinfo:phytree:plot:UnknownTypeName',...
                                'Unknown option for %s.',upper(okargs{k}));
                        else
                            if l==4
                                l=1;
                            end
                            renderType = oktypes{l};
                        end
                    case 2 % orientation
                        oktypes={'left','right','top','bottom'};
                        l = strmatch(lower(pval),oktypes); %#ok
                        if isempty(l)
                            error('Bioinfo:phytree:plot:UnknownOrientation',...
                                'Unknown option for %s.',upper(okargs{k}));
                        else
                            orientation = oktypes{l};
                        end
                    case 3 % rotation
                        if isreal(pval(1))
                            rotation = double(pval(1));
                        else
                            error('Bioinfo:phytree:plot:NotValidType',...
                                'ROTATION must be numeric and real');
                        end
                    case 4 % the leaf order func 'LeafOrder'
                        LeafOrder = pval;
                    case 5 % colors of the leafs 'cellColor'
                        cellColorsVector = pval;
                    case 6 % colors of the branches 'clusterColor'
                        clusterColorsVector = pval;
                    case 7 % size of the branches'clusterSize'
                        clusterSizeVector = pval;
                    case 8 % 'intermediateNodeSize'
                        intermediateNodeSize = pval;
                    case 9 % 'intermediatoNodeLabel'
                        intermediatoNodeLabel = pval;
                    case 10 % labels of the leafs 'LeafLabels'
                        leafLabels = pval;
                    case 11 % labels of the branches 'branchLabels'
                        branchLabels = pval;
                    case 12 % display leaf labels
                        if pval,
                            display_leaf_labels=tree;
                        else
                            display_leaf_labels=false;
                        end
                    case 13 % (Noa)'Group'
                        group = pval;   
                    case 14 % display branch labels
                        if pval,
                            display_branch_labels='on';
                        else
                            display_branch_labels='off';
                        end
                end
            end
        end
    end
    %%
    % set dependent defaults
    if isnan(display_branch_labels)
        if isequal(renderType,'radial')
            display_branch_labels = true;
        else
            display_branch_labels = false;
        end
    end
    if isnan(display_leaf_labels)
        if isequal(renderType,'radial')
            display_leaf_labels = true;
        else
            display_leaf_labels = false;
        end
    end

    %%

    tree = doBasicCalculations(tree,activeBranches,renderType);
    x=tree.x;
    y=tree.y;
    nodeIndex   = 1:tree.numLabels;
    leafIndex   = 1:tree.numLeaves;
    branchIndex = tree.numLeaves+1:tree.numLabels;

    %%
    % check empty names
    for ind = nodeIndex
        if isempty(tree.names{ind})
            if ind > tree.numLeaves
                %tr.names{ind} = ['Branch ' num2str(ind-tr.numLeaves)];
            else
                tree.names{ind} = ['Leaf ' num2str(ind)];
            end
        end
    end

    %%
    % rendering graphic objects
    fig = figure('Renderer','ZBuffer');
    h.axes = axes; hold on;
    sepUnit = max(tree.x)*[-1/20 21/20];


    %%
    % setting the axes
    switch renderType
        case {'square','angular'}
            switch orientation
                case 'left'
                    set(h.axes,'YTick',1:numel(tree.terminalNodes),'Ydir','reverse',...
                        'YtickLabel','','YAxisLocation','Right')
                    if display_leaf_labels
                        set(h.axes,'Position',[.05 .10 .7 .85])
                    else
                        set(h.axes,'Position',[.05 .10 .9 .85])
                    end
                    xlim(sepUnit);
                    ylim([0 numel(tree.terminalNodes)+1]);
                case 'right'
                    set(h.axes,'YTick',1:numel(tree.terminalNodes),'Xdir','reverse','Ydir','reverse',...
                        'YtickLabel','','YAxisLocation','Left')
                    if display_leaf_labels
                        set(h.axes,'Position',[.25 .10 .7 .85])
                    else
                        set(h.axes,'Position',[.05 .10 .9 .85])
                    end
                    xlim(sepUnit);
                    ylim([0 numel(tree.terminalNodes)+1]);
                case 'top'
                    set(h.axes,'XTick',1:numel(tree.terminalNodes),...
                        'XtickLabel','','XAxisLocation','Top')
                    if display_leaf_labels
                        set(h.axes,'Position',[.10 .05 .85 .7])
                    else
                        set(h.axes,'Position',[.10 .05 .85 .9])
                    end
                    ylim(sepUnit);
                    xlim([0 numel(tree.terminalNodes)+1]);
                case 'bottom'
                    set(h.axes,'XTick',1:numel(tree.terminalNodes),'Xdir','reverse','Ydir','reverse',...
                        'XtickLabel','','XAxisLocation','Bottom')
                    if display_leaf_labels
                        set(h.axes,'Position',[.10 .25 .85 .7])
                    else
                        set(h.axes,'Position',[.10 .05 .85 .9])
                    end
                    ylim(sepUnit);
                    xlim([0 numel(tree.terminalNodes)+1]);
            end
        case 'radial'
            set(h.axes,'XTick',[],'YTick',[])
            set(h.axes,'Position',[.05 .05 .9 .9])
            display_leaf_labels = false;
            axis equal
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % drawing lines
    switch renderType
        case 'square'
            X = tree.x([nodeIndex;repmat([tree.par(1:tree.numLabels-1) tree.numLabels],2,1)]);
            Y = tree.y([repmat(nodeIndex,2,1);[tree.par(1:tree.numLabels-1) tree.numLabels]]);
            switch orientation
                case {'left','right'}
                    h.BranchLines = plot(X,Y,'-k');
                    delete(h.BranchLines(~tree.activeNodes))
                    h.BranchLines = h.BranchLines(tree.activeNodes);
                case {'top','bottom'}

                    if scale==LOGSCALE,
                        Y(3,end-1) = Y(2,end-1);
                        X(2,end-1) = X(1,end-1)/2;
                        X(3,end-1) = X(1,end-1)/2;
                        plot(Y(:,2:end),X(:,2:end),'-k');
                    else                

                    h.BranchLines = plot(Y,X,'-k');
                    delete(h.BranchLines(~tree.activeNodes))
                    h.BranchLines = h.BranchLines(tree.activeNodes);
                    end
            end
        case 'angular'
            X = tree.x([nodeIndex;[tree.par(1:tree.numLabels-1) tree.numLabels]]);
            Y = tree.y([nodeIndex;[tree.par(1:tree.numLabels-1) tree.numLabels]]);
            switch orientation
                case {'left','right'}
                    h.BranchLines = plot(X,Y,'-k');
                    delete(h.BranchLines(~tree.activeNodes))
                    h.BranchLines = h.BranchLines(tree.activeNodes);
                case {'top','bottom'}
                    h.BranchLines = plot(Y,X,'-k');
                    delete(h.BranchLines(~tree.activeNodes))
                    h.BranchLines = h.BranchLines(tree.activeNodes);
            end
        case 'radial'
            R = tree.x;
            A = tree.y / numel(tree.terminalNodes)*2*pi+rotation*pi/180;
            tree.x = R .* sin(A);
            tree.y = R .* cos(A);
            X = tree.x([nodeIndex;[tree.par(1:tree.numLabels-1) tree.numLabels]]);
            Y = tree.y([nodeIndex;[tree.par(1:tree.numLabels-1) tree.numLabels]]);
            h.BranchLines = plot(X,Y,'-k');
            delete(h.BranchLines(~tree.activeNodes))
            h.BranchLines = h.BranchLines(tree.activeNodes);
    end

    %%
    % drawing default nodes (with no coloring)
    switch renderType
        case {'square','angular'}
            switch orientation
                case {'left','right'}
                    h.LeafDots = plot(tree.x(leafIndex(tr.activeNodes(leafIndex))),...
                        tree.y(leafIndex(tree.activeNodes(leafIndex))),'o',...
                        'MarkerSize',MARKER_SIZE,'MarkerEdgeColor','k',...
                        'MarkerFaceColor','k');
                case {'top','bottom'}
                    if scale==LOGSCALE,
                        Y_to_plot = tree.y(leafIndex(tree.activeNodes(leafIndex)));
                        Y_to_plot(1) = Y(1,end-1);
                        X_to_plot = tree.x(leafIndex(tree.activeNodes(leafIndex)));
                        X_to_plot(1) = X(2,end-1);
                        h.LeafDots = plot(Y_to_plot,X_to_plot,'o','MarkerSize',MARKER_SIZE,'MarkerEdgeColor','k',...
                            'MarkerFaceColor','k');
                    else

                        h.LeafDots = plot(tree.y(leafIndex(tree.activeNodes(leafIndex))),...
                            tree.x(leafIndex(tree.activeNodes(leafIndex))),'o',...
                            'MarkerSize',MARKER_SIZE,'MarkerEdgeColor','k',...
                            'MarkerFaceColor','k');
                    end

            end
        case 'radial'
            h.BranchDots = plot(tr.x(branchIndex(tr.activeNodes(branchIndex))),...
                tr.y(branchIndex(tr.activeNodes(branchIndex))),'o',...
                'MarkerSize',MARKER_SIZE,'MarkerEdgeColor','k',...
                'MarkerFaceColor','k');
            h.LeafDots = plot(tr.x(leafIndex(tr.activeNodes(leafIndex))),...
                tr.y(leafIndex(tr.activeNodes(leafIndex))),'square',...
                'MarkerSize',MARKER_SIZE,'MarkerEdgeColor','k',...
                'MarkerFaceColor','w');
    end

    %%
    % remove underscores from groups

    for iter=1:length(group)-2,
        for ind = 1:length(group{iter})
            group{iter}{ind}(group{iter}{ind}=='_')=' ';
        end
    end
    
    %%
    % color group nodes appropriately

    for iter=1:length(group)-2,
        terminal_color=terminal_colors{iter};
        for i=1:length(tree.names),
            ind=find(strcmp(group{iter},tree.names{i}));
            if ~isempty(ind),
                if length(terminal_color)==3,
                    if strcmp(orientation,'bottom')
                    plot(tree.y(i),tree.x(i),['b' 'o'],'MarkerSize',MARKER_SIZE,'MarkerEdgeColor',terminal_color,...
                        'MarkerFaceColor',terminal_color);
                else
                    plot(tree.x(i),tree.y(i),['b' 'o'],'MarkerSize',MARKER_SIZE,'MarkerEdgeColor',terminal_color,...
                        'MarkerFaceColor',terminal_color);
                    end
                else
                if strcmp(orientation,'bottom')
                    plot(tree.y(i),tree.x(i),[terminal_color 'o'],'MarkerSize',MARKER_SIZE,'MarkerEdgeColor',terminal_color,...
                        'MarkerFaceColor',terminal_color);
                else
                    plot(tree.x(i),tree.y(i),[terminal_color 'o'],'MarkerSize',MARKER_SIZE,'MarkerEdgeColor',terminal_color,...
                        'MarkerFaceColor',terminal_color);
                end
                end
            end
        end
    end

    % % color significant branch nodes in green
    % if ~isempty(group{end})
    % sig_branch_size = zeros(1,length(group{end}));
    % sig_branch_size(find(group{end}<0.00001)) = 5;
    % sig_branch_size(find(group{end}>=0.00001 & group{end}<0.0001)) = 4;
    % sig_branch_size(find(group{end}>=0.0001 & group{end}<0.001)) = 3;
    % sig_branch_size(find(group{end}>=0.001)) = 2;
    % for i=1:length(tr_names2),
    %     ind=find(strcmp(group{end-1},tr_names2{i}));
    %     if ~isempty(ind),
    %         if strcmp(orientation,'bottom')
    %              %plot(tr.y(i),tr.x(i),[internal_color 'o'],'MarkerSize',MARKER_SIZE,'MarkerEdgeColor',internal_color,...
    %              %    'MarkerFaceColor',internal_color);  
    %              plot(Y(:,i),X(:,i),['-' internal_color], 'LineWidth', sig_branch_size(ind));             
    %         else
    %               %plot(tr.x(i),tr.y(i),[internal_color 'o'],'MarkerSize',MARKER_SIZE,'MarkerEdgeColor',internal_color,...
    %               %   'MarkerFaceColor',internal_color); 
    %              plot(X(:,i),Y(:,i),['-' internal_color], 'LineWidth', sig_branch_size(ind));
    %              
    %              
    %              
    %         end
    %     end
    % end
    % end

    indicator_names = zeros(1,length(tr_names2));
    if ~isempty(group{end})
        groups = group{end};
        groups_names = group{end-1};
        for k=1:length(groups),
            cur_group = groups{k};
            cur_group_names = groups_names{k};
            sig_branch_size = zeros(1,length(cur_group));
            sig_branch_size(find(cur_group<0.00001)) = 5;
            sig_branch_size(find(cur_group>=0.00001 & cur_group<0.0001)) = 4;
            sig_branch_size(find(cur_group>=0.0001 & cur_group<0.001)) = 3;
            sig_branch_size(find(cur_group>=0.001)) = 2;
            for i=1:length(tr_names2),
                ind=find(strcmp(cur_group_names,tr_names2{i}));
                if ~isempty(ind),
                    if indicator_names(i)==0,
                        line_type = '-';
                    elseif indicator_names(i)==1,
                        line_type = '--';
                    elseif indicator_names(i)==2,
                        line_type = ':';
                    end
                    indicator_names(i) = indicator_names(i)+1;
                    if strcmp(orientation,'bottom')
                        %plot(Y(:,i),X(:,i),['-' internal_color{k}], 'LineWidth', sig_branch_size(ind));
                        plot(Y(:,i),X(:,i),line_type,'color',internal_color{k}, 'LineWidth', sig_branch_size(ind));
                    else
                        %plot(X(:,i),Y(:,i),['-' internal_color{k}], 'LineWidth', sig_branch_size(ind));
                        plot(X(:,i),Y(:,i),line_type,'color', internal_color{k}, 'LineWidth', sig_branch_size(ind));

                    end
                    hold on;
                end
            end
        end
    end

    %%
    % resize figure if needed
    switch renderType
        case {'square','angular'}
            switch orientation
                case {'left','right'}
                    correctFigureSize(fig, 15 * numel(tr.terminalNodes),0);
                    fontRatio = max(get(fig,'Position').*[0 0 0 1])/numel(tr.terminalNodes);
                case {'top','bottom'}
                    correctFigureSize(fig, 0, 15 * numel(tr.terminalNodes));
                    fontRatio = max(get(fig,'Position').*[0 0 1 0])/numel(tr.terminalNodes);
            end
        case 'radial'
            temp = 10/pi*numel(tr.terminalNodes);
            correctFigureSize(fig,temp,temp);
            fontRatio = max(get(fig,'Position').*[0 0 1 0])/numel(tr.terminalNodes);
    end

    set(h.axes,'Fontsize',min(9,ceil(fontRatio/1.5)));


    % set branch node labels
    if dispBranchLabels,
        X = tr.x(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels)));
        Y = tr.y(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels)));  
        switch renderType
            case {'square','angular'}
                switch orientation
                    case {'left'}
                        if scale==LINSCALE,
                            h.branchNodeLabels = text(X+sepUnit(1)/2,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                         elseif scale==LOGSCALE,
                            h.branchNodeLabels = text(X,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                        end
                        set(h.branchNodeLabels,'color',[0 0.75 0.75],'clipping','on')
                        set(h.branchNodeLabels,'vertical','bottom')
                        set(h.branchNodeLabels,'horizontal','right')
                       %set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                       set(h.branchNodeLabels,'Fontsize',7);
                    case {'right'}
                        h.branchNodeLabels = text(X+sepUnit(1)/2,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                        set(h.branchNodeLabels,'color',[0 0.75 0.75],'clipping','on')
                        set(h.branchNodeLabels,'vertical','bottom')
                        set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                    case {'top'}
                        h.branchNodeLabels = text(Y,X-sepUnit(1)/2,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                        set(h.branchNodeLabels,'color',[0 0.75 0.75],'clipping','on')
                        set(h.branchNodeLabels,'vertical','bottom','Rotation',30)
                        set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                    case {'bottom'}
                        if scale==LINSCALE,
                            h.branchNodeLabels = text(Y,X+sepUnit(1)/2,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                        elseif scale==LOGSCALE,
                            h.branchNodeLabels = text(Y+0.5,X,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                        end
                        %set(h.branchNodeLabels,'color',[0 0.7 1],'clipping','on')
                        set(h.branchNodeLabels,'color','k','clipping','on')
                        set(h.branchNodeLabels,'vertical','bottom','Rotation',0)
                        %set(h.branchNodeLabels,'vertical','bottom')
                        set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                        set(h.branchNodeLabels,'Fontsize',8);
                end
            case 'radial'
                h.branchNodeLabels = text(X,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
                set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
                set(h.branchNodeLabels,'vertical','bottom')
                set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio*1.2)));
                for ind = 1:numel(h.branchNodeLabels)
                    if X(ind)<0
                        set(h.branchNodeLabels(ind),'horizontal','right')
                        set(h.branchNodeLabels(ind),'Position',get(h.branchNodeLabels(ind),'Position')+[sepUnit(1)/2 0 0])
                    else
                        set(h.branchNodeLabels(ind),'horizontal','left')
                        set(h.branchNodeLabels(ind),'Position',get(h.branchNodeLabels(ind),'Position')-[sepUnit(1)/2 0 0])
                    end
                end
        end
    end
    % set leaf nodes labels
    if dispLeafLabels,
        X = tr.x(leafIndex(tr.activeNodes(1:tr.numLeaves)));
        Y = tr.y(leafIndex(tr.activeNodes(1:tr.numLeaves)));
        switch renderType
            case {'square','angular'}
                switch orientation
                    case {'left'}
                        h.leafNodeLabels = text(X-sepUnit(1)/2,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                        %set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                        set(h.leafNodeLabels,'color','b','clipping','on')
                        set(h.leafNodeLabels,'horizontal','left')
                        set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                    case {'right'}
                        h.leafNodeLabels = text(X-sepUnit(1)/2,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                        set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                        set(h.leafNodeLabels,'horizontal','right')
                        set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                    case {'top'}
                        h.leafNodeLabels = text(Y,X-sepUnit(1)/2,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                        set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                        set(h.leafNodeLabels,'horizontal','left','Rotation',60)
                        set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                    case {'bottom'}
                        h.leafNodeLabels = text(Y,X-sepUnit(1),tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                        set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                        set(h.leafNodeLabels,'horizontal','right','Rotation',60)
                        set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));
                end
            case 'radial'
                h.leafNodeLabels = text(X,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
                set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
                set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio*1.2)));
                % textHeight = mean(cell2mat(get(h.leafNodeLabels,'Extent')))*[0 0 0 1]';
                for ind = 1:numel(h.leafNodeLabels)
                    if X(ind)<0
                        set(h.leafNodeLabels(ind),'horizontal','right')
                        set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')+[sepUnit(1) 0 0])
                    else
                        set(h.leafNodeLabels(ind),'horizontal','left')
                        set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')-[sepUnit(1) 0 0])
                    end
                    %             a=atan(Y(ind)/X(ind))*180/pi;
                    %             if a > 0  a = max(0,a-60)/2; else
                    %                       a = min(0,a+60)/2; end
                    %             set(h.leafNodeLabels(ind),'Rotation',a)
                end
                [sortedY,hsY]=sort(Y);
                idx=hsY(X(hsY)>0 & sortedY>0);
                if numel(idx)
                    extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
                    positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
                    for i = 2:numel(idx)
                        position = get(h.leafNodeLabels(idx(i)),'Position');
                        positionY = max(positionY+extentY,position(2));
                        position(2) = positionY;
                        set(h.leafNodeLabels(idx(i)),'Position',position)
                    end
                end
                idx=hsY(X(hsY)<0 & sortedY>0);
                if numel(idx)
                    extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
                    positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
                    for i = 2:numel(idx)
                        position = get(h.leafNodeLabels(idx(i)),'Position');
                        positionY = max(positionY+extentY,position(2));
                        position(2) = positionY;
                        set(h.leafNodeLabels(idx(i)),'Position',position)
                    end
                end
                idx=flipud(hsY(X(hsY)>0 & sortedY<0));
                if numel(idx)
                    extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
                    positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
                    for i = 2:numel(idx)
                        position = get(h.leafNodeLabels(idx(i)),'Position');
                        positionY = min(positionY-extentY,position(2));
                        position(2) = positionY;
                        set(h.leafNodeLabels(idx(i)),'Position',position)
                    end
                end
                idx=flipud(hsY(X(hsY)<0 & sortedY<0));
                if numel(idx)
                    extentY = get(h.leafNodeLabels(idx(1)),'Extent')*[0;0;0;1];
                    positionY = get(h.leafNodeLabels(idx(1)),'Position')*[0;1;0];
                    for i = 2:numel(idx)
                        position = get(h.leafNodeLabels(idx(i)),'Position');
                        positionY = min(positionY-extentY,position(2));
                        position(2) = positionY;
                        set(h.leafNodeLabels(idx(i)),'Position',position)
                    end
                end
        end
    end


    % correct axis limits given the extent of labels
    if dispBranchLabels
        E = cell2mat(get(h.branchNodeLabels,'Extent'));
        if strcmp(get(gca,'XDir'),'reverse')
            E(:,1) = E(:,1) - E(:,3);
        end
        if strcmp(get(gca,'YDir'),'reverse')
            E(:,2) = E(:,2) - E(:,4);
        end
        E=[E;[xlim*[1;0] ylim*[1;0] diff(xlim) diff(ylim)]];
        mins = min(E(:,[1 2]));
        maxs = max([sum(E(:,[1 3]),2) sum(E(:,[2 4]),2)]);
        axis([mins(1) maxs(1) mins(2) maxs(2)])
    end




    if dispLeafLabels
        E = cell2mat(get(h.leafNodeLabels,'Extent'));
        if strcmp(get(gca,'XDir'),'reverse')
            E(:,1) = E(:,1) - E(:,3);
        end
        if strcmp(get(gca,'YDir'),'reverse')
            E(:,2) = E(:,2) - E(:,4);
        end
        E=[E;[xlim*[1;0] ylim*[1;0] diff(xlim) diff(ylim)]];
        mins = min(E(:,[1 2]));
        maxs = max([sum(E(:,[1 3]),2) sum(E(:,[2 4]),2)]);
        axis([mins(1) maxs(1) mins(2) maxs(2)])
    end

    % set terminal nodes labels
    switch renderType
        case {'square','angular'}
            X = tr.x(tr.terminalNodes) * 0;
            Y = tr.y(tr.terminalNodes);
            switch orientation
                case {'left'}
                    X = X + max(xlim) - sepUnit(1)/2;
                    h.terminalNodeLabels = text(X,Y,tr.names(tr.terminalNodes),'Color','b');
                case {'right'}
                    X = X + max(xlim) - sepUnit(1)/2;
                    h.terminalNodeLabels = text(X,Y,tr.names(tr.terminalNodes));
                    set(h.terminalNodeLabels,'Horizontal','right')
                case {'top'}
                    X = X + max(ylim) - sepUnit(1)/2;
                    h.terminalNodeLabels = text(Y,X,tr.names(tr.terminalNodes));
                    set(h.terminalNodeLabels,'Rotation',90)
                case {'bottom'}
                    X = X + max(ylim) - sepUnit(1)/2;
                    h.terminalNodeLabels = text(Y,X,tr.names(tr.terminalNodes));
                    set(h.terminalNodeLabels,'Rotation',270)
            end
        case 'radial'
            h.terminalNodeLabels = text(0,0,' ');
    end

    % set group text color
    if dispTerminalLabels
        for iter=1:length(group)-2,
            terminal_color=terminal_colors{iter};
            for i=1:length(tr.names(tr.terminalNodes)),
                ind=find(strcmp(group{iter},tr.names{tr.terminalNodes(i)}));
                if ~isempty(ind),
                    set(h.terminalNodeLabels(i),'color',terminal_color,'clipping','off')
                end
            end
        end
    end

    if scale==LOGSCALE,
        set(h.axes,'Yscale','log');
        set(gca,'YLimMode','auto');
    end

    set(gca,'XTick',[]);



    % if dispTerminalLabels
    %     set(h.terminalNodeLabels,'Fontsize',min(9,ceil(fontRatio/1.5)));
    % end

    % if ~dispBranchLabels
    %     set(h.branchNodeLabels,'visible','off');
    % end
    % if ~dispLeafLabels
    %     set(h.leafNodeLabels,'visible','off');
    % end
    if ~dispTerminalLabels
        set(h.terminalNodeLabels,'visible','off');
    end

    box on
    hold off

    % store handles
    set(fig,'UserData',h)
    if nargout
        handles = h;
    end

    %set(h.branchNodeLabels,'visible',display_branch_labels);


