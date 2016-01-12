function entropy = Tree_Entropy(Tree,groups)

[Leaves_codes] = prepare_leaves_codes_for_entropy_score(Tree,groups);

poi=get(Tree,'Pointers');
n=get(Tree,'NUMLEAVES');

p=0.99;
un=unique(Leaves_codes);
if(un(1)==-1)
    un=un(2:end);
end
for l1=1:length(un)
    for l2=l1+1:length(un)
        type1=un(l1);
        type2=un(l2);        
        
        poiA=[zeros(n,2);poi];
        poiA=[[1:2*n-1]' poiA];        
        
        for i=2*n-1:-1:n+1%Writing the dircet children of each node
            poiA(poiA(i,2),4)=i;
            poiA(poiA(i,3),4)=i;
        end
        
        for i=1:n%writing the anscestors of each leaf
            anst=[];
            noot=0;
            loc=i;
            while noot==0
                anst=[anst poiA(loc,4)];
                loc=poiA(loc,4);
                if loc==2*n-1
                    noot=1;
                end
                ansc{i}=anst;
            end
        end        
        
        child=cell(2*n-1,1);
        for i=1:n%Writing the descendant of each node
            ind=1;
            noot=0;
            while noot==0
                loc=ansc{i}(ind);
                child{loc}=[child{loc} i];
                ind=ind+1;
                if loc==2*n-1
                    noot=1;
                end
            end
        end
        
        poiA(find(Leaves_codes==type1),5)=1;
        poiA(find(Leaves_codes==type2),6)=1;
        for i=n+1:2*n-1%Running on the inner nodes and writing the dcendant of eahc of the two types each node has
            an=child{i};
            n1=length(find(Leaves_codes(an)==type1));
            n2=length(find(Leaves_codes(an)==type2));
            poiA(i,5)=n1;
            poiA(i,6)=n2;
        end
        
        for i=1:2*n-1
            if poiA(i,5)>p*(poiA(i,5)+poiA(i,6))
                poiA(i,7)=type1;
                poiA(i,8)=poiA(i,5)+poiA(i,6);
            end
            if poiA(i,6)>p*(poiA(i,5)+poiA(i,6))
                poiA(i,7)=type2;
                poiA(i,8)=poiA(i,5)+poiA(i,6);
            end
            
        end
        
        poiA(:,9:10)=poiA(:,7:8);
        for i=1:n
            an=ansc{i};
            fd=find(poiA(an,7)~=0);
            if length(fd)>0
                poiA(i,9:10)=0;
            end
            if length(fd)>1
                poiA(an(fd(1:end-1)),9:10)=0;
            end
        end
        
        clus1=length(find(poiA(:,9)==type1));
        clus2=length(find(poiA(:,9)==type2));
        n1=length(find(Leaves_codes==type1));
        n2=length(find(Leaves_codes==type2));
        
        permu = log(nchoosek(max(n1-clus1,clus1-1),min(n1-clus1,clus1-1))) + log(factorial(clus1)) +...
            log(nchoosek(max(n2-clus2,clus2-1),min(n2-clus2,clus2-1))) + log(factorial(clus2));
        %prob=permu/(factorial(n1)*factorial(n2));
        entropy(l1,l2)=permu;entropy(l2,l1)=entropy(l1,l2);
    end
end
end


function [Leaves_codes] = prepare_leaves_codes_for_entropy_score(T,groups)

leaves_names = get(T,'LeafNames');
Leaves_codes = repmat(-1,1,length(leaves_names));
for i = 1:length(groups),
    group  = groups{i};
    for j=1:length(group),
        ind = find(strcmp(leaves_names,group{j})==1);
        Leaves_codes(ind) = i;
    end
end

end



