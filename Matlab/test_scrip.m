addpath('S:\Noa\Tree_analysis_Feb2015\Ver_3\');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Step 1: Filer Data %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MS_tab_input_file = 'S:\Noa\Tree_analysis_Feb2015\Projects\Halaban\YUCLAT\MS_table.tab'; % Samples x Loci   YUCLAT
SNP_tab_input_file = 'S:\Noa\Tree_analysis_Feb2015\Projects\Halaban\YUCLAT\snp_table.txt'; % Samples x SNPs
MS_tab_output_file = 'IO_Files\MS_table_filtered';
SNP_tab_output_file = 'IO_Files\SNP_table_filtered.txt';
PERCENT_LOCI_AMPLIFIED_PER_CELL = 0.3;
PERCENT_CELLS_AMPLIFIED_PER_LOCUS = 0.3;

Filter_Data('MStableInputFile',MS_tab_input_file,'MStableOutputFile',MS_tab_output_file,...
    'SNPtableInputFile',SNP_tab_input_file,'SNPtableOutputFile',SNP_tab_output_file,...
    'PercentLociAmplifiedPerCell',PERCENT_LOCI_AMPLIFIED_PER_CELL,'PercentCellsAmplifiedPerLocus',PERCENT_CELLS_AMPLIFIED_PER_LOCUS);


%Filter_Data('MStableInputFile','/net/mraid11/export/data/dcsoft/home/Noa/Tree_analysis_Feb2015/Projects/Halaban/YUCLAT/MS_table.tab','MStableOutputFile','/net/mraid11/export/data/dcsoft/home/Noa/Tree_analysis_Feb2015/Ver_3/IO_Files/MS_output_tab.txt')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Step 2: Root calculation %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MS_tab_input_file = MS_tab_output_file;
SNP_tab_input_file = SNP_tab_output_file;
MS_root_output_file = 'IO_Files\MS_root_identifier.txt';
SNP_root_output_file = 'IO_Files\SNP_root_identifier.txt';
ROOT = 'AVE'; %can be 'Ave' or list of cell names

Root_Calculation('MStableInputFile',MS_tab_input_file,'SNPtableInputFile',SNP_tab_input_file,...
    'MSrootOutputFile',MS_root_output_file,'SNProotOutputFile',SNP_root_output_file,'RootCalculation',ROOT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Step 3: Distance MAtrix %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MS_root_input_file = MS_root_output_file;
SNP_root_input_file = SNP_root_output_file;
Distance_matrix_output_file = 'IO_Files\Distance_matrix.txt';
METRIC = 0;
CELLS_TO_ANALYZED = 'MS_SNP'; %MS/SNP/MS_SNP/list of cells to be analyzed
DATA_TO_RECONSTRUCT = 'MS_SNP';
MS_WEIGHT = 0.001;

Distance_Calculation('MStableInputFile',MS_tab_input_file,'SNPtableInputFile',SNP_tab_input_file,...
    'MSrootInputFile',MS_root_input_file,'SNProotInputFile',SNP_root_input_file,'DistanceMatOutputFile',Distance_matrix_output_file,...
   'DistMetric',METRIC,'UseDataToReconstructTree',DATA_TO_RECONSTRUCT,'CellsToBeAnalysed',CELLS_TO_ANALYZED,'MSweight',MS_WEIGHT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Step 4: Reconstructing algorithm %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Distance_Matrix_input_file = Distance_matrix_output_file;
Tree_output_file_newick = 'IO_Files\Tree_output.newick';
Tree_output_file_matlab = 'IO_Files\Tree_output.tree';
Rec_Algo = 'NJ';

Tree_Distance_Based_Reconstruction('DistanceMatInputFile',Distance_Matrix_input_file,'Algorithm',Rec_Algo,...
    'TreeOutputFileNewick',Tree_output_file_newick,'TreeOutputFileMatlab',Tree_output_file_matlab);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Step 5: Tree Bootstraping %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootIterations_val = 1;
Bootstrap_output_file = 'IO_Files\Bootstrap_vals.txt';
Tree_input_file_newick = Tree_output_file_newick;

Tree_Bootstraping('MStableInputFile',MS_tab_input_file,'SNPtableInputFile',SNP_tab_input_file,...
    'MSrootInputFile',MS_root_input_file,'SNProotInputFile',SNP_root_input_file,'BootIterations',BootIterations_val,'BootstrapValsOutputFile',Bootstrap_output_file,...
   'DistMetric',METRIC,'UseDataToReconstructTree',DATA_TO_RECONSTRUCT,'CellsToBeAnalysed',CELLS_TO_ANALYZED,'MSweight',MS_WEIGHT,...
   'Algorithm',Rec_Algo,'TreeInputFileNewick',Tree_input_file_newick);



Bootstrap_input_file = Bootstrap_output_file;
Newick_output_file = 'IO_Files\fuzzy_NJ_newick.txt';
Fuzzy_NJ('TreeInputFileNewick',Tree_input_file_newick,'BootstrapValsInputFile',Bootstrap_input_file,'NewickOutputFile',Newick_output_file);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Step 5: Tree Scoring %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell_data_input_file = 'C:\Users\noac\Desktop\Halaban_cell_data.csv';
Tree_input_file_newick = Tree_output_file_newick;



%Tree_Scoring('CellDataInputFile',cell_data_input_file,'TreeInputFileNewick',Tree_input_file_newick,);




disp('end');