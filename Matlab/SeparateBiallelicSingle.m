function [MSTableSep, MSNamesSep, LeafNames] = SeparateBiallelicSingle(CallFile1,MaxMutNum)

Data1 = importdata(CallFile1);

LeafNames = Data1.textdata(2:end,1);
MSNames = Data1.textdata(1,2:end);
MS1 = Data1.data;

c=1;
MSTableSep = [];
MSNamesSep = {};

for i=1:size(MS1,2)
    ms = MS1(:,i);
    len = length(ms);
    Table = [ms];
    if (sum(~isnan(Table)) < 2)
        continue;
    end
       
    [a b] = count_unique(Table);
    S = zeros(50,1);
    S(a) = b;
    [pck locs] = findpeaks(S);
    [Idx,D] = knnsearch(locs,Table);
    Table(D > MaxMutNum) = nan;
    for k=1:length(locs)
        tab = nan(size(Table));
        tab(Idx == k) = Table(Idx == k);
        tab1 = tab(1:len);
        inds = isnan(tab1);
        A=sum(~isnan(tab1)>0);
        
        if (sum(~isnan(tab1))<2)
            continue;
        end
        MSTableSep = [MSTableSep tab1];
        MSNamesSep{c} = [MSNames{i} '_' num2str(k)];
        c = c+1;
    end
end
