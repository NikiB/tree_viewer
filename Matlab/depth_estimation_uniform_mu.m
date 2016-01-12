function t=depth_estimation_uniform_mu(ML_table,val)

% 6/1/07
% slices the relevant vector from ML and returnes its maximum t (naive).
% For several loci it computes the likelihood of each depth and returnes
% the maximum, the mean, std and the actual likelihood vector (normalized
% for probabilities).

% filter out loci where there is no information on mutation rate/bias/value



val_range=-((size(ML_table,1)-1)/2):(size(ML_table,1)-1)/2;
t_range=1:size(ML_table,2);

logP=zeros(1,size(ML_table,2));
for i=1:length(val),
    if ~isnan(val(i)),
        val_ind=find(val_range==val(i));
        if val(i)>max(val_range),
            val_ind=length(val_range);
        end
        if val(i)<min(val_range),
            val_ind=1;
        end
        logP=logP+log(ML_table(val_ind,:));
    end
end
[y,ii]=max(logP);
t=t_range(ii);