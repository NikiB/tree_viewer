function [D,num_common_loci,ind_common_loci]=calculate_matrix_distances(mat,metric,weights)

% Input:
%       - mat: mutation matrix (cells X loci)
%       - metric:
%               0 for Absolute
%               1 for Normalized Absolute
%               2 for Euclidean
%               3 for Squared
%               4 for Mutation count (Equal or Not)
%               5 for Maximum Likelihood estimator
%
% Output:
%       - D: distance mateix (cells X cells)
%       - num_common_loci: cells X cells matrix which contains the number of common loci that were compared for each pair of cells.
%       - ind_common_loci: cells X cells matrix which contains the list of
%       the common loci.

if nargin<3,
    weights = ones(1,size(mat,2));
end

if metric==5,
    [ML_tab] = create_ML_table(1/100);
    %[ML_tab] = create_ML_table([ - 3.49979857e-03, 3.63701136e-04]);
end



NAN_IDENTIFIER=-100;




D = zeros(size(mat,1));
num_common_loci = zeros(size(mat,1));
ind_common_loci = {};
for i=1:size(mat,1),
    for j=(i+1):size(mat,1),
        indin=find(~isnan(mat(i,:))&(~isnan(mat(j,:))));
        num_common_loci(i,j)=length(indin);
        ind_common_loci{i,j} = indin;
        cur_weights = weights(indin);
        if isempty(indin),
            D(i,j) = NaN;
        else
            if metric==0,
                % Absolute distance
                A = mat(i,:)-mat(j,:);
                B = abs(A);
                D(i,j) = nanmean(B);
                
                d_temp= mean(abs(round((mat(i,indin)-mat(j,indin)))));
                if D(i,j)~=d_temp,
                    disp('D not equal');
                end
              
            elseif metric==1
                % Normalized Absolute distance
                A = mat(i,indin)/max(1,sum(abs(mat(i,indin))));
                B = mat(j,indin)/max(1,sum(abs(mat(j,indin))));
                D(i,j) = mean(abs(A-B));
            elseif metric==2,
                % Euclidean distance
                D(i,j) = sqrt(sum((round((mat(i,indin)-mat(j,indin)))).^2));
            elseif metric==3
                % Squared distance
                D(i,j)= mean((round((mat(i,indin)-mat(j,indin)))).^2);
            elseif metric==4,
                % Mutation count (Equal or Not) distance
                inds=find(round(mat(i,indin)-mat(j,indin))~=0);
                vec = zeros(1,length(indin));
                vec(inds)=1;  
                vec = vec.*cur_weights;
                D(i,j) = sum(vec) / length(indin);
                %D(i,j)=length(find(round(mat(i,indin)-mat(j,indin))~=0))/length(indin);
            elseif metric==5,
                %Maximum Likelihood estimator
                
                vec=round((mat(i,indin)-mat(j,indin)));
                D(i,j)=depth_estimation_uniform_mu(ML_tab,vec);
            end
            
            
        end
    end 
end

D = D+D';
num_common_loci = num_common_loci+num_common_loci';
for k=2:length(ind_common_loci)
    for p=1:k-1
        ind_common_loci{k,p} = ind_common_loci{p,k};
    end
end

% finally when D is NaN assign the average distance
d=[];
for i=1:size(D,1),
    d = [d D(i,i+1:end)];
end
D(isnan(D)) = mean(d(~isnan(d)));

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ML_tab] = create_ML_table(mu)

X = -15:15; 
T = 1:1000;

 
%Single step (only +1 or -1)
xstep=-3:3;
middle=(length(xstep)-1)/2+1;
pstep=zeros(1,length(xstep));
pstep(middle+1)=0.5;
pstep(middle-1)=0.5;

ML_tab = build_ML_table(X,T,mu,xstep,pstep);
save ML_tab ML_tab;

end


                
              