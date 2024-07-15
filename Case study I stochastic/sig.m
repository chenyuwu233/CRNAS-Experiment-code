% This is a function that calculate the gradient of variance function sig such
% that sig(p,beta,nu,b,E,m) = p (beta + nu - H)/(beta - nu + H)(exp(2(beta - nu + H) t) - exp((beta - nu + H)t))
%
% Input:
%   - d: (scalar) dosage levels
%   - t: (scalar) time points 
%   - theta: (1 x 6): (p,beta,nu,b,E,m)
%   - cmd: (string) command that what information required: 'grad', 'Hess'
%
% Output:
%   - obj: (scalar)
%   - grad: (NT x ND x 6) [p,beta,nu,b,E,m]
%   - Hess: (6 x 6 matrix) [p,beta,nu,b,E,m]
%
% Requirment:
%   - mu(d,t,theta)
%   - ratio(d,beta,nu,b,E,m)

function [obj,grad,Hess] = sig(d,t,theta,cmd)
    %% Assign
    p = theta(1);
    beta = theta(2);
    nu = theta(3);
    b = theta(4);
    E = theta(5);
    m = theta(6);

    %%
    theta_np = theta;
    theta_np(1) = 1;
    [obj_e_np, grad_e_np, Hess_e_np] = mu(d,t,theta_np,cmd);
    obj_E = p.*(obj_e_np^2 - obj_e_np);
    grad_e_np = squeeze(grad_e_np);
    Hess_e_np = squeeze(Hess_e_np);
    grad_e_np(1) = [];
    Hess_e_np(1,:) = [];
    Hess_e_np(:,1) = [];
    grad_E = zeros(1,6);
    grad_E(2:6) = p.*(2*obj_e_np-1).*grad_e_np;
    grad_E(1) = (obj_e_np^2 - obj_e_np);
    Hess_E = zeros(6,6);
    Hess_E(2:6,2:6) = p.*(2.*(grad_e_np)*(grad_e_np') + (2*obj_e_np-1).*Hess_e_np);
    Hess_E(1,2:6) = (2*obj_e_np-1).*grad_e_np;
    Hess_E(2:6,1) = (2*obj_e_np-1).*grad_e_np;

 
    [obj_r,grad_r,Hess_r] = ratio(d,beta,nu,b,E,m);




    
    %% Objective
    obj = obj_r*obj_E;

    %% Gradient
    grad = zeros(1,6);
    grad(1) = obj_r.*grad_E(1);
    grad(2:6) = obj_E.*grad_r + obj_r.*grad_E(2:6);


    %% Hessian
    Hess = zeros(6,6);
    Hess(2:6,2:6) = obj_E.*Hess_r + (grad_E(2:6)')*grad_r + grad_r'*grad_E(2:6) + obj_r.*Hess_E(2:6,2:6);
    Hess(1,2:6) = grad_E(1).*grad_r + obj_r.*Hess_E(1,2:6);
    Hess(2:6,1) = Hess(1,2:6); 
    
end