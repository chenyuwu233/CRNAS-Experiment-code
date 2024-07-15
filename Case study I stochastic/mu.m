% This is a function that calculate the gradient of mean function mu such
% that mu(p,beta,nu,b,E,m) = p exp((beta-nu+H(d;b,E,m))t)
%
% Input:
%   - d: (1 x ND) dosage levels
%   - t: (1 x NT) time points 
%   - theta: (1 x 6): (p,beta,nu,b,E,m)
%   - cmd: (string) command that what information required: 'grad', 'Hess'
%
% Output:
%   - obj: (scalar)
%   - grad: (NT x ND x 6) [p,beta,nu,b,E,m]
%   - Hess: (6 x 6 matrix) [p,beta,nu,b,E,m]
%
% Requirment:
%   - grad_pop(vec,d,t,cmd)


function [obj,grad,Hess] = mu(d,t,theta,cmd)
    %% Assign
    p = theta(1);
    beta = theta(2);
    nu = theta(3);
    b = theta(4);
    E = theta(5);
    m = theta(6);

%%    
    vec1 = [beta - nu,b,E,m];
    if cmd == 'Hess'
        [Pop,grad_lam,grad_b,grad_E,grad_m,Hess_temp] = grad_pop(vec1,d,t,1);
        obj  = p.*Pop;
        
        grad = cat(3,Pop,p.*grad_lam,-p.*grad_lam,p.*grad_b,p.*grad_E,p.*grad_m);
        Hess = zeros(length(t),length(d),6,6);
        Hess(:,:,4:6,4:6) = p.*Hess_temp(:,:,2:4,2:4);
        Hess(:,:,2,2) = p.*Hess_temp(:,:,1,1); 
        Hess(:,:,3,3) = p.*Hess_temp(:,:,1,1); 
        Hess(:,:,2,3) = -p.*Hess_temp(:,:,1,1);
        Hess(:,:,3,2) = -p.*Hess_temp(:,:,1,1);
        Hess(:,:,2,4:6) = p.*Hess_temp(:,:,1,2:4); % Birth
        Hess(:,:,4:6,2) = p.*Hess_temp(:,:,2:4,1); % Birth
        Hess(:,:,3,4:6) = -p.*Hess_temp(:,:,1,2:4); % Death
        Hess(:,:,4:6,3) = -p.*Hess_temp(:,:,2:4,1); % Death
        Hess(:,:,1,1) = 0;
        Hess(:,:,1,2:6) = cat(3,grad_lam,-grad_lam,grad_b,grad_E,grad_m);
        Hess(:,:,2:6,1) = cat(3,grad_lam,-grad_lam,grad_b,grad_E,grad_m);
    else
        [Pop,grad_lam,grad_b,grad_E,grad_m] = grad_pop(vec1,d,t);
        obj  = p.*Pop;
        grad = cat(3,Pop,p.*grad_lam,-p.*grad_lam,p.*grad_b,p.*grad_E,p.*grad_m);
        Hess = zeros(6,6);
    end
end