function bootstraps=bootstrap_vals(partition_real,partitions_bootstrap)

% bootstrap_vals=bootstrap_vals(partition_real,partitions_bootstrap)
% 21/1/07
% for each node (partition) returnes the bootstrap value
% Input:
%       partition_real - an N*1 vector of length N nodes (including
%                        internal nodes)
%       partitions_bootstrap - an N*m matrix where each column is of length N nodes (including
%                        internal nodes) and m is the number of bootstrap
%                        iteration

m=size(partitions_bootstrap,2);
bootstraps=zeros(length(partition_real),1);
part_all=partitions_bootstrap(:);
for i=1:length(partition_real),
    ind=find(part_all==partition_real(i));
    if isempty(ind),
        L=0;
    else
        L=length(ind);
    end
    bootstraps(i)=L;
end
bootstraps=floor(100*bootstraps/m);
bootstraps(bootstraps>100)=100;
