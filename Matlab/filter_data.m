function [MS_mat,SNP_mat,cell_names,cell_indin] = filter_data(MS_mat,SNP_mat,cell_names,params)
    
if ~isempty(MS_mat),
    [MS_indin,cells_MS_indin] = filter_mat(MS_mat,params);
end
if ~isempty(SNP_mat),
    [SNP_indin,cells_SNP_indin] = filter_mat(SNP_mat,params);
end

if strcmp(params.DATA_TO_RECONSTRUCT,'MS_SNP'),
    cell_indin = intersect(cells_MS_indin,cells_SNP_indin);
    MS_mat = MS_mat(cell_indin,MS_indin);
    SNP_mat = SNP_mat(cell_indin,SNP_indin);
    cell_names = cell_names(cell_indin);
elseif strcmp(params.DATA_TO_RECONSTRUCT,'MS'),
    MS_mat = MS_mat(cells_MS_indin,MS_indin);
    cell_names = cell_names(cells_MS_indin);
    cell_indin = cells_MS_indin;
elseif strcmp(params.DATA_TO_RECONSTRUCT,'SNP'),
    SNP_mat = SNP_mat(cells_SNP_indin,SNP_indin);
    cell_names = cell_names(cells_SNP_indin);
    cell_indin = cells_SNP_indin;
end
  


end


function [loci_indin,cells_indin] = filter_mat(mat,params)
sm = sum(isnan(mat'));
if params.PERCENT_LOCI_AMPLIFIED_PER_CELL>1,
   smnn = sum(~isnan(mat'));
   cells_indin = find(smnn>=params.PERCENT_LOCI_AMPLIFIED_PER_CELL); 
else
    cells_indin = find(sm/size(mat,2)<=(1-params.PERCENT_LOCI_AMPLIFIED_PER_CELL));
end

sm = sum(isnan(mat));
loci_indin=find(sm/size(mat,1)<=(1-params.PERCENT_CELLS_AMPLIFIED_PER_LOCUS));

end
