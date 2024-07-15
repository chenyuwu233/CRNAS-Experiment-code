%% gradient for pL/(1+exp(-cv*n+nx0))

function grad = logistic_grad_sub(p,L,n,nx0,cv)
    grad = zeros(4,1);
    exp_temp = exp(-cv*n+nx0);
    grad(1) = L/(1+exp_temp);
    grad(2) = p/(1+exp_temp);
    grad(3) = L*p*cv*exp_temp/(1+exp_temp)^2;
    grad(4) = -L*p*exp_temp/(1+exp_temp)^2;
end