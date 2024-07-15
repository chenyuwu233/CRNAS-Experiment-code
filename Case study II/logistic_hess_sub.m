%% Hessian for pL/(1+exp(-cv*n+nx0))

function hess = logistic_hess_sub(p,L,n,nx0,cv)
    hess = zeros(4,4);
    exp_temp = exp(-cv*n+nx0);
    hess(1,2) = 1/(1+exp_temp);
    hess(1,3) = L*cv*exp_temp/(1+exp_temp)^2;
    hess(1,4) = -L*exp_temp/(1+exp_temp)^2;
    hess(2,3) = p*cv*exp_temp/(1+exp_temp)^2;
    hess(2,4) = -p*exp_temp/(1+exp_temp)^2;
    hess(3,3) = p*L*cv*(-cv*exp_temp*(1+exp_temp) + 2*cv*exp_temp^2)/(1+exp_temp)^3;
    hess(3,4) = -p*L*cv*(exp_temp*(1+exp_temp) - 2*exp_temp)/(1+exp_temp)^3;
    hess(4,4) = -p*L*(exp_temp*(1+exp_temp) - 2*exp_temp)/(1+exp_temp)^3;
    hess(2,1) = hess(1,2);
    hess(3,1) = hess(1,3);
    hess(4,1) = hess(1,4);
    hess(3,2) = hess(2,3);
    hess(4,2) = hess(2,4);
    hess(4,3) = hess(3,4);
end