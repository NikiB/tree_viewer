function [pval,pval_red,pval_white]=hypergeometric_statistics(N,red,white,s,xred,xwhite)

pval_red=sum(hygepdf(xred:min(s,red),N,s,red));
pval_white=sum(hygepdf(xwhite:min(s,white),N,s,white));
if xred/s>=red/N,
    pval=pval_red;
end
if xred/s<red/N,
    pval=pval_white;    
end
pval=pval_red;

    

