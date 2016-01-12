% MAIN SCRIPT FOR TREE RECONSTRUCTION AND SNP ANALYSIS

addpath('S:\Noa\Tree_analysis_Feb2015\Ver_2.1\');

% MS_file = 'S:\LINEAGE\Hiseq\NSR4\fastq_human\Calling\NSR4_AC_X_mat_1a__0_01__30.tab'; 
% cell_data_file = 'S:\Noa\Rambam_Collaboration\NSR4\';
% SNP_file = [];
% 
% %first we need to take only the relevant cells from the MS file
% [num char raw] = xlsread(cell_data_file);
% data.cellID = num;
% dataMS = importdata(SNP_file);
% cell_names_MS = dataMS.textdata(2:end,1);
% MS_names = dataMS.textdata(1,2:end);
% MS_mat = dataMS.data;
% cell_indin = [];
% for i=1:length(cell_names_MS),
%     cur_cell_name = cell_names_MS(i);
%     indID = findstr(cur_cell_name,'cID');
%     curr_ID = str2num(cur_cell_name(indID+1:indID+5));
%     if find(data.cellID==curr_ID),
%         cell_indin = [cell_indin; i];
%     end
% end
% 
% cID2393


MS_file = 'S:\Noa\Tree_analysis_Feb2015\Projects\Halaban\YUCLAT\MS_table.tab'; % Samples x Loci   YUCLAT
SNP_file = 'S:\Noa\Tree_analysis_Feb2015\Projects\Halaban\SNP_table_full.txt'; % Samples x SNPs
cell_data_file = 'S:\Noa\Tree_analysis_Feb2015\Projects\Halaban\YUCLAT\Halaban_cell_data.csv';
% cell_data_file = 'C:\Users\noac\Desktop\Halaban_cell_data.csv';

Results_dir = 'S:\Noa\Tree_analysis_Feb2015\Projects\Halaban\YUCLAT\Results\';

TreeParams.ALGORITHM = 'NJ'; 
TreeParams.METRIC = 4; % 0-Absolute/1-Normalized Absolute/2-Euclidean/3-Squared/4-Mutation count (Equal or Not)
TreeParams.PERCENT_LOCI_AMPLIFIED_PER_CELL = 0.1; % values from 0 to 1
TreeParams.PERCENT_CELLS_AMPLIFIED_PER_LOCUS = 0.3; % values from 0 to 1
TreeParams.MS_WEIGHT = 0.01;
TreeParams.DATA_TO_RECONSTRUCT = 'MS'; % MS/SNP/MS_SNP
TreeParams.CELLS_TO_ANALYZED = 'MS'; % MS/SNP/MS_SNP
TreeParams.ROOT = 'AVE';
TreeParams.INFORMATIVE_CELL_NAMES_CRITERION = [];%'Sample Name'; % can be empty or one of the titles in the cell data file
TreeParams.BOOT_ITER = 1;
TreeParams.VALIDATE_ITER = 0;



Lineage_Analysis(MS_file,SNP_file,cell_data_file,Results_dir,TreeParams);
