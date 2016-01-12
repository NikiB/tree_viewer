function [] = Lineage_Analysis(MS_file,SNP_file,Cell_data_file,Results_dir,TreeParams,figure_name)

 
% Read the MS mutation calling file
MS_mat = [];
if ~isempty(MS_file),
    %[MS_mat, MS_names, cell_names_MS] = SeparateBiallelicSingle(MS_file,2);
    dataMS = my_importdata(MS_file);
    cell_names_MS = dataMS.textdata(2:end,1);
    MS_names = dataMS.textdata(1,2:end);
    MS_mat = dataMS.data;
end

% Read the SNP mutation calling file
SNP_mat = [];
if ~isempty(SNP_file),
    dataSNP = importdata(SNP_file);
    cell_names_SNP = dataSNP.textdata(2:end,1);
    SNP_names = dataSNP.textdata(1,2:end);
    SNP_mat = dataSNP.data;
end

fig_name = [Results_dir figure_name  '_' TreeParams.DATA_TO_RECONSTRUCT '_' num2str(TreeParams.METRIC) '_'...
  num2str(TreeParams.PERCENT_LOCI_AMPLIFIED_PER_CELL*100) '_' num2str(TreeParams.PERCENT_CELLS_AMPLIFIED_PER_LOCUS*100)];  

%intersect cell names 
if strcmp(TreeParams.CELLS_TO_ANALYZED,'MS_SNP'),    
   [cell_names,ind_MS,ind_SNP] = intersect(cell_names_MS,cell_names_SNP);
   MS_mat = MS_mat(ind_MS,:);
   SNP_mat = SNP_mat(ind_SNP,:);  
elseif strcmp(TreeParams.CELLS_TO_ANALYZED,'MS'),
    cell_names = cell_names_MS;        
elseif strcmp(TreeParams.CELLS_TO_ANALYZED,'SNP'),
    cell_names = cell_names_SNP;    
end


if strcmp(TreeParams.DATA_TO_RECONSTRUCT,'MS_SNP'),    
   %add the MS weight to the fig name 
   s = num2str(TreeParams.MS_WEIGHT);
   ind_point = findstr(s,'.');
   w = ['0' s(ind_point+1:end)];
   fig_name = [fig_name '_MSweight_' w]; 
end

%read the cell_data file
cell_data = read_cell_data_file(Cell_data_file);

%reconstruct tree
[Tree] = Tree_Analysis(MS_mat,SNP_mat,cell_names,cell_data,TreeParams,fig_name);



if ~isempty(SNP_file),
    cell_names_SNP = generate_short_names(cell_names_SNP,cell_data,TreeParams);
    SNP_Analysis(Tree,SNP_mat,cell_names_SNP,SNP_names,fig_name);
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = read_cell_data_file(cell_data_file_name)

% [pathstr, name, ext]= fileparts(cell_data_file_name);
% 
% Data = importdata(cell_data_file_name);
% 
% if strcmp(ext,'.csv'),
%    data.data_types = strsplit(Data{1},','); 
% else
%     data.data_types = strsplit(Data{1},'\t');
% end
% 
% for i=1:length(data.data_types),
%     data.all_data{i} = [];
% end
% 
% for i=2:length(Data),
%     if strcmp(ext,'.csv'),
%         info = strsplit(Data{i},',');
%     else
%         info = strsplit(Data{i},'\t');
%     end
%     for j=1:length(info),
%         data.all_data{j} =  [data.all_data{j} info(j)];
%     end
% end

[num char raw] = xlsread(cell_data_file_name);
data.cellID = num;

data.data_types = raw(1,:);

for i=1:length(data.data_types),
     data.all_data{i} = [];
end

for i=1:length(data.data_types),
     data.all_data{i} = raw(2:end,i);
end

end

