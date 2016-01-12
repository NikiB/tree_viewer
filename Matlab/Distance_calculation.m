function [D] = Distance_calculation(MS_mat,SNP_mat,params)

% Input:
%       - MS_mat:   MS mutation matrix (cells X loci)
%       - SNP_mat:  SNP mutation matrix (cells X SNPs)
%       - params: tree reconstruction parameter file
% Output:
%       - D: distance matrix (cells X cells)


data_to_reconstruct = params.DATA_TO_RECONSTRUCT; % MS/SNP/MS_SNP



%if (params.MS_FLAG & params.SNP_FLAG)
if strcmp(data_to_reconstruct,'MS_SNP'),
    unite_mat = [MS_mat SNP_mat];
    weights = [params.MS_WEIGHT*ones(1,size(MS_mat,2)) ones(1,size(SNP_mat,2))];
    [D] = calculate_matrix_distances(unite_mat,4,weights);
    return;
end


%if params.MS_FLAG
if strcmp(data_to_reconstruct,'MS'),
    weights = ones(1,size(MS_mat,2));
   [D_MS] = calculate_matrix_distances(MS_mat,params.METRIC,weights);
end

%if params.SNP_FLAG
if strcmp(data_to_reconstruct,'SNP'),
    weights = ones(1,size(SNP_mat,2));
    [D_SNP] = calculate_matrix_distances(SNP_mat,4,weights);
end

D = [];
% if (params.MS_FLAG & params.SNP_FLAG)
%     D = params.MS_WEIGHT*D_MS + D_SNP; 
% elseif params.MS_FLAG,
%     D = D_MS;
% elseif params.SNP_FLAG,
%     D = D_SNP;
% end

if strcmp(data_to_reconstruct,'MS_SNP'),
    %D = params.MS_WEIGHT*D_MS + D_SNP; 
    D = D_MS + D_SNP; 
elseif strcmp(data_to_reconstruct,'MS'),
    D = D_MS;
elseif strcmp(data_to_reconstruct,'SNP'),
    D = D_SNP;
end




%   -add ML calculation

end





