function ML_table=build_ML_table(x,T,mu,xstep,pstep)

% ML_table=build_ML_table(x,T,mu,xstep,pstep)
% Builds a probability table where ML_table(x,T) is the probability to see
% a mutation of size x after T divisions, given that in a single step the
% probability of a mutation happening is mu, and if it happens the
% probability of it being xstep is pstep.


MAX_LEN = 500; % Max length of MS

if (length(mu) == 1)
    ML_table=zeros(length(x),length(T));
    
    for j=1:length(x),
        if x(j)==0,
            ML_table(j,1)=1-mu;
        else
            ind=find(xstep==x(j));
            if isempty(ind),
                ML_table(j,1)=0;
            else
                ML_table(j,1)=mu*pstep(ind);
            end
        end
    end


    for i=2:size(ML_table,2),
        for jnew=1:length(x),
            add=0;
            for jold=1:length(x),
                deltax=x(jnew)-x(jold);
                ind=find(xstep==deltax);
                if ~isempty(ind),
                    add=add+pstep(ind)*ML_table(jold,i-1);
                end
            end
            ML_table(jnew,i)=(1-mu)*ML_table(jnew,i-1)+mu*add;
        end
    end
end

% mu is dependent on 2 parameters
if (length(mu) == 2)
    P = zeros(MAX_LEN);
    len = length(P);
    Mu0 = mu(1);
    Mu1 = mu(2);
    Mu = Mu0 + Mu1*(1:len);
    
    for j=1:length(xstep),
        if xstep(j)==0,
            P(1:len+1:len*len) = 1-Mu;
        else
            if (xstep(j) > 0)
                StartInd = xstep(j)+1;
                P(StartInd:len+1:len*(len-xstep(j))) = Mu(abs(xstep(j))+1:len)*pstep(j);
            else
                StartInd = len*abs(xstep(j))+1;
                P(StartInd:len+1:len*len-xstep(j)) = Mu(1:len-abs(xstep(j)))*pstep(j);
            end
            %P(2:len+1:len*len) = Mu(2:len)/2;
            %P(len+1:len+1:len*len) = Mu(1:len-1)/2;
        end
    end
    
    ML_table = cell(T,1);
    ML_table{1} = P;

    for i=2:length(T)
        ML_table{i} = ML_table{i-1}*P;
    end    
end