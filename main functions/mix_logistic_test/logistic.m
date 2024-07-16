function ret = logistic(x,L,n,nx0)
    
    ret = L/(1+exp(-n*x+nx0));
end