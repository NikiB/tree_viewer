function [t,Z] = seqneighjoin(D, method, names, varargin) 
%SEQNEIGHJOIN neighbor-joining method for phylogenetic tree reconstruction.
%
%  TREE = SEQNEIGHJOIN(DIST) computes a phylogenetic tree object from the
%  pairwise distances DIST between the species or products applying the
%  neighbor-joining method. The input DIST is a matrix (or vector) such as
%  is generated by SEQPDIST.
%
%  TREE = SEQNEIGHJOIN(DIST,METHOD) selects the method to compute the
%  distances of the new nodes to all other nodes at every iteration. The
%  general expression to calculate the distances between the new node (n)
%  (after joining i and j) and all others nodes (k) is given by: 
%
%      D(n,k) =  a*D(i,k) + (1-a)*D(j,k) - a*D(n,i) - (1-a)*D(n,j) 
%
%  This expression guarantees to find the correct tree with additive data
%  (a.k.a. minimum variance reduction). The options for METHOD are:
%
%     'equivar'    --- Assumes equal variance and independence of
%     (default)        evolutionary distance estimates (a = 1/2). Such as
%                      in Studier and Keppler, JMBE (1988). 
%     'firstorder' --- Assumes a first order model of the variances and
%                      covariances of evolutionary distance estimates, 'a'
%                      is adjusted at every iteration to a value between 0
%                      and 1. Such as in Gascuel, JMBE (1997). 
%     'average'    --- New distances are the weighted average of previous
%                      distances, the branch distances are ignored:
%                            D(n,k) =  [ D(i,k) + D(j,k) ] /2
%                      As in the original neighbor-joining algorithm by
%                      Saitou and Nei, JMBE (1987)
%
%  TREE = SEQNEIGHJOIN(DIST,METHOD,NAMES) passes a list of names to label
%  the leaf nodes (e.g. species or products) in the phylogenetic tree
%  object. NAMES can be a vector of structures with the fields 'Header' or
%  'Name' or a cell array of strings. In both cases the number of elements
%  provided must comply with the number of samples used to generate the
%  pairwise distances in DIST. 
%
%  TREE = SEQNEIGHJOIN(...,'REROOT',false) excludes rerooting the resulting
%  tree. This is useful to observe the original linkage order followed by
%  the algorithm. By default SEQNEIGHJOIN reroots the resulting tree using
%  the mid-point method. 
%
% Example:
%
%     % Load a multiple alignment of amino acids:
%     seqs = fastaread('pf00002.fa');
% 
%     % Measure the 'Jukes-Cantor' pairwise distances:
%     dist = seqpdist(seqs,'method','jukes-cantor','indels','pair');
% 
%     % Build the phylogenetic using the neighbor-joining algorithm
%     tree = seqneighjoin(dist,'equivar',seqs)
%     view(tree)
%
% See also MULTIALIGN, PHYTREE, PHYTREE/REROOT, PHYTREE/VIEW, SEQLINKAGE,
% SEQPDIST.  

% References: 
% [1] Saitou N, Nei M.The neighbor-joining method: a new method for
%     reconstructing phylogenetic trees. Mol Biol Evol.(1987) 4(4):406-25
% [2] Gascuel O. BIONJ: An improved version of the NJ algorithm based on a
%     simple model of sequence data. Mol. Biol. Evol. (1997) 14:685-695  
% [3] Studier JA, Keppler KJ. A note on the neighbor-joining algorithm of 
%     Saitou and Nei. Mol Biol Evol. (1988) 5(6):729-31. 

% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.1.8.1 $  $Date: 2005/06/09 21:57:44 $

rerootTree = false;
checkForNames = true;

% check the input distances
if isnumeric(D)
    [m, n] = size(D);
    if isvector(D)
        n = numel(D);
        m = (1+sqrt(1+8*n))/2; % number of leaves
        if m ~= fix(m)
            error('Bioinfo:seqneighjoin:DbadSize',...
                'Size of DIST not compatible with the output of the SEQPDIST function.');
        end
        D = squareform(D);
    elseif m~=n
        error('Bioinfo:seqneighjoin:DnotSquare',...
              'Size of DIST not compatible with the output of the SEQPDIST function.');
    end
else
    error('Bioinfo:seqneighjoin:DnotNumeric',...
          'DIST must be a numeric vector compatible with the output of the SEQPDIST function.');
end

% Selects appropiate method
if nargin == 1   
    method = 'e'; % set default method
    checkForNames = false;
else
    okmethods = {'equivar','firstorder','average','reroot'};
    methodkeys = {'e','f','a'};  
    k = find(strncmpi(method,okmethods,numel(method)));
    if isempty(k)
        error('Bioinfo:seqneighjoin:UnknownMethod',...
              'Unknown method name: %s.',method);
    elseif length(k)>1
        error('Bioinfo:seqneighjoin:IncorrectMethod',...
              'Ambiguous method name: %s.',method);
    elseif k<4
        method = methodkeys{k};
    else % case that second and third input are optional paired inputs
        if numel(varargin) || nargin==2
            error('Bioinfo:seqneighjoin:IncorrectInputArguments',...
                  'Incorrect format of input arguments.')
        end
        varargin = {method,names};
        checkForNames = false;
        method = 'e'; % set default method
    end
end

% detects the names
if checkForNames && (nargin >2) % names were supplied, check validity
    if iscell(names) || isfield(names,'Header') || isfield(names,'Name') || isfield(names,'LocusName')
        if isfield(names,'Header')  % if struct put them in a cell
            names = {names(:).Header};
        elseif isfield(names,'Name') % if struct put them in a cell
            names = {names(:).Name};
        elseif isfield(names,'LocusName') % if struct put them in a cell
            names = {names(:).LocusName};
        end
        names = names(:);
        namesSupplied = true;
        if numel(names)~=m
            error('Bioinfo:seqneighjoin:IncorrectSize',...
            'NAMES must have the same size as number of leaves in the tree')
        end
    elseif strncmpi(names,{'reroot'},numel(names))
        varargin = {names varargin{:}};
        namesSupplied = false;        
    else
        error('Bioinfo:seqneighjoin:IncorrectInputType',...
          'NAMES must be a cell with char arrays or a vector of structures.')
    end
else
    namesSupplied = false;
end

% check optional input arguments
 nvarargin = numel(varargin);
 if nvarargin 
     rerootTree = true;
     root_name=varargin{2};
 end
%     if rem(nvarargin,2)
%         error('Bioinfo:seqneighjoin:IncorrectNumberOfArguments',...
%               'Incorrect number of arguments to %s.',mfilename);
%     end
%     okargs = {'reroot',''};
%     for j=1:2:nvarargin
%         pname = varargin{j};
%         pval = varargin{j+1};
%         k = find(strncmpi(pname,okargs,numel(pname)));
%         if isempty(k)
%             error('Bioinfo:seqneighjoin:UnknownParameterName',...
%                 'Unknown parameter name: %s.',pname);
%         elseif length(k)>1
%             error('Bioinfo:seqneighjoin:AmbiguousParameterName',...
%                 'Ambiguous parameter name: %s.',pname);
%         else
%             rerootTree = opttf(pval);
%             if isempty(rerootTree)
%                 error('Bioinfo:seqneighjoin:rerootInputOptionNotLogical',...
%                        '%s must be a logical value, true or false.',...
%                        upper(char(okargs(k))));
%             end
%         end
%     end
% end

% ------------------------------------------------------------------------
% Algorithm starts here:
%
%   D       -  pairwise distance matrix (full form)
%   method  -  'a' for Saitou & Nei, 'e' for Studier & Keppler, and 'f' for
%              Gascuel
%   names   -  labels for leaf nodes

N  = size(D,1);
Z  = zeros(N-1,2);
bd = zeros(N*2-1,1); 

p = 1:N;   % pointers in matrix
bc=1;      % branch counter
if method=='f'
    V=D;   %initialize the variance matrix for Gascuel method
end
    
for n = N:-1:3
    R = sum(D)/(n-2);                 % sums of columns 
    Q = D-repmat(R,n,1)-repmat(R',1,n)+diag(inf(n,1)); %Studier & Keppler optimized
    % S = (sum(sum(D))/(n-2)+Q)/2     % original total sum in Saitou & Nei  
    [m,g] = min(Q(:));  
    [i,j] = ind2sub(n,g);             % find minimum
    if i>j k=i;i=j;j=k; end           % j>i always
    pp = p([i j]);                    % pointers to join
    bl = (R([i j])*[1 -1;-1 1] + D(j,i))/2; % branch lengths
    bl = max(bl,0);
    bd(pp) = bl;                      % save branch lengths
    Z(bc,:) = pp;                     % save pointers
    h = [1:i-1 i+1:j-1 j+1:n];
    switch method                     % distances to new node
        case 'a' % Saitou & Nei method
           d = (sum(D(h,[i j]),2))/2;      
           D = [[D(h,h) d];[d' 0]];   % update distance matrix
        case 'e' % Studier & Keppler method
           d = (sum(D(h,[i j]),2)-D(j,i))/2; 
           if any(d<0)
               d = max(d,0);
           end
           D = [[D(h,h) d];[d' 0]];   % update distance matrix
        case 'f' % Gascuel method
           if V(i,j)
              lambda = max(0,min(1,(1+sum(V(h,i)-V(h,j))/(n-2)/V(i,j))/2));
           else
              lambda = 1/2;
           end
           d = D(h,[i j])*[lambda;1-lambda] - bl*[lambda;1-lambda];
           if any(d<0)
               d = max(d,0);
           end
           D = [[D(h,h) d];[d' 0]];   % update distance matrix
           v = V(h,[i j])*[lambda;1-lambda] - V(i,j)*lambda*(1-lambda);
           if any(v<0)
               v = max(v,0);
           end
           V = [[V(h,h) v];[v' 0]];   % update variance matrix
    end
    p  = [p(h),N+bc];                 % update pointers
    bc = bc+1;                        % update branch counter
end
Z(bc,:) = p;                          % pointers of last branch
bd(p) = D(2)/2;                       % lengths of last branch

% convert data to a phylogenetic tree object
if namesSupplied
    t = phytree(Z,bd,names);
else
    t = phytree(Z,bd);
end
if rerootTree
    temp=get(t,'LeafNames');
    root_ind=find(strcmp(temp,root_name));
    t = reroot(t,root_ind,0); % reroot with mid-point method
end