function [T,groups] = Tree_Analysis(T,MS_mat,cell_names,cell_data)
    
    Tree_to_Tsne(T,MS_mat,cell_names);
    groups = Groups_Analysis(T,cell_data);

    [T,groups] = Enrichment_Calculation(T,groups);
    

end