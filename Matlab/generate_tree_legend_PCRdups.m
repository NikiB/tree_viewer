function [  ] = generate_tree_legend_PCRdups(Params,fig_name)

fid = fopen([fig_name '.txt'],'w');
fprintf(fid,'Tree reconstructed algorithm: %s\n',Params.ALGORITHM);
fprintf(fid,'Distace matrix: '); 
switch Params.METRIC
    case 0
        fprintf(fid,'Absolute');
    case 1
        fprintf(fid,'Normalized Absolute');
    case 2
        fprintf(fid,'Euclidean');
    case 3
        fprintf(fid,'Squared');
    case 4
        fprintf(fid,'Mutation count (Equal or Not)');
    case 5
        fprintf(fid,'Maximum Likelihood estimator');
end
fprintf(fid,'\n');
fprintf(fid,'Minimum percent of amplified loci for each cell: %d\n', Params.PERCENT_LOCI_AMPLIFIED_PER_CELL*100);
fprintf(fid,'Minimum percent of cells which were amplified for each locus: %d\n', Params.PERCENT_CELLS_AMPLIFIED_PER_LOCUS*100);
fprintf(fid,'Root calculation: ');

if strcmp(Params.ROOT,'AVE')
       fprintf(fid,'average signature of all cells\n');
else
   fprintf(fid,'average signature of the following cell types: ');
   for i=1:length(Params.ROOT)-1,
      fprintf(fid,'%s, ',Params.ROOT{i}); 
   end
   fprintf(fid,'%s\n',Params.ROOT{end}); 
end


fprintf(fid,'Colored leaves represent PCR duplicates, each color for each pair of duplicates\n');
fclose all;

end

