function [T] = tree_validation_test(T,D,MS_mat,SNP_mat,cell_names,root_name,params)


%validate tree
leaves_amount = get(T,'NumLeaves');
Tr_validate = [];
Dr_validate = [];
MS_matr_validate = [];
SNP_matr_validate = [];

iter_amount = max(floor(leaves_amount*0.50),10);
for i=1:iter_amount,
    %remove i leaves randomally
    for j=1:params.VALIDATE_ITER,
        [vals, inds] = sort(rand(1,leaves_amount-1)); % it is -1 because of the root
        indin = setdiff([1:leaves_amount],inds(1:i));
        MS_matr = MS_mat(indin,:);
        SNP_matr = SNP_mat(indin,:);
        
        Dr = Distance_calculation(MS_matr,SNP_matr,params);
    
        %reconstruct the tree
        if isempty(params.ROOT),
            [Tr,Z] = my_seqneighjoin(Dr,'equivar',cell_names);
        else
            [Tr,Z] = my_seqneighjoin(Dr,'equivar',cell_names,'reroot',root_name);
        end
            
        Tr_validate{i,j} = Tr;
        Dr_validate{i,j} = Dr;
        MS_matr_validate{i,j} = MS_matr;
        SNP_matr_validate{i,j} = SNP_matr;
    end    
end

end
