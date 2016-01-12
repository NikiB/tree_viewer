function [pval,alpha_val,N,s1,s2,x]=tree_enrichments(t,group_names,q)


% [pval,alpha_val,N,s1,s2,x]=tree_enrichment(t,group_names)
% Given a tree and names of leafs belonging to a group, test for each
% branch the hypergeometric enrichment of the group within the branch
% children.
% alpha_val is determined according to FDR of q

if nargin<3,
    q=0.2;
end

leaf_names=get(t,'LeafNames');

[children,ancs]=get_branch_children(t);

N=length(leaf_names);
s1=length(group_names);
s2=zeros(1,length(children));
x=zeros(1,length(children));
pval=zeros(1,length(children));

for i=1:length(children),
    set2test=children{i};
    s2(i)=length(set2test);
    x(i)=length(intersect(set2test,group_names));
    % remove root
    if (s2(i)~=N)
        pval(i)=hypergeometric_statistics(N-1,s1,N-1-s1,s2(i),x(i),s2(i)-x(i));
    else
        pval(i)=2;
    end
end
alpha_val=fdr(pval,q);



