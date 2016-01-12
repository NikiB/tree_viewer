% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = doBasicCalculations(tr,activeBranches,renderType)

    % helper function to compute and find some features of the tree
    tr.numLeaves = tr.numBranches + 1;
    tr.numLabels = tr.numBranches + tr.numLeaves; 


    % remove uderscores from names
    for ind = 1:tr.numLabels
        tr.names{ind}(tr.names{ind}=='_')=' ';
    end

    % obtain parents for every node
    tr.par(tr.tree(:)) = tr.numLeaves + [1:tr.numBranches 1:tr.numBranches];

    % find active nodes
    tr.activeNodes = true(tr.numLabels,1);
    for ind =tr.numBranches:-1:1
        tr.activeNodes(tr.tree(ind,:)) = tr.activeNodes(tr.numLeaves+ind) & activeBranches(ind);
    end

    % propagate last leaf
    tr.lastleaf = 1:tr.numLabels;
    for ind = tr.numBranches:-1:1
        if ~tr.activeNodes(tr.tree(ind,1))
            tr.lastleaf(tr.tree(ind,:))=tr.lastleaf(ind+tr.numLeaves);
        end
    end

    tr.activeBranches = tr.activeNodes(tr.numLeaves+1:tr.numLabels)&activeBranches;
    tr.activeLeaves = tr.activeNodes(1:tr.numLeaves);

    % find x coordinates of branches
    tr.x = tr.dist;
    for ind = tr.numBranches:-1:1
        tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
    end

    % find y coordinates of branches
    tr.terminalNodes = tr.lastleaf([true,diff(tr.lastleaf(1:tr.numLeaves))~=0]);
    tr.y=zeros(tr.numLabels,1);
    tr.y(tr.terminalNodes)=1:length(tr.terminalNodes);
    switch renderType
        case 'square'
            for ind = 1:tr.numBranches
                if tr.activeBranches(ind)
                    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
                end
            end
        case {'angular','radial'}
            for ind = 1:tr.numBranches
                if tr.activeBranches(ind)
                    if tr.x(tr.tree(ind,1))/tr.x(tr.tree(ind,2))>3
                        tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,1));
                    elseif tr.x(tr.tree(ind,2))/tr.x(tr.tree(ind,1))>3
                        tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,2));
                    else
                        tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
                    end
                end
            end
    end
end