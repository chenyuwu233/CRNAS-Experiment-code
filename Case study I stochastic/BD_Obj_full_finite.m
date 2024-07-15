% This is an objective function with Branching process, for 2 subpopulation
%
% 
% Input: 
%   - Params: (1 x 13 vector) all necessary parameters ([p_i,r_i, d_i, b_i, E_i^n_i, n_i], c)
%   - DATA: (NT x ND x NR matrix) observation
%   - Conc: (1 x ND vector) concentration levels
%   - Time: (1 x NT vector) of time points
%   - n: (scalar) initial total cell number.
%   - eps_grad: (scalar) epsilon for the gradient finite difference
%   - eps_hess: (scalar) epsilon for the hessian finite difference
%
% Output:
%   - obj: (scalar) objective value
%   - grad: (1 x 6*S+1 vector) gradient at Params
%   - Hess: (6*S+1 x 6*S+1 matrix) Hessian matrix at Params
%
% Requirment:
%   - mu(d,t,theta)
%   - sig(d,t,theta)

function [obj,grad,Hess] = BD_Obj_full_finite(Params,DATA,Conc,Time,n,eps_grad,eps_hess)    
    obj = BD_Obj_full(Params,DATA,Conc,Time,n,'Obj');
    grad = zeros(1,length(Params));
    Hess = zeros(length(Params),length(Params));
    vec  = zeros(size(Params));
    for i = 1:length(Params)
        ei = vec;
        ei(i) = 1;
        % ol = BD_Obj_full(Params - eps_grad*ei,DATA,Conc,Time,n,'Obj');
        or = BD_Obj_full(Params + eps_grad*ei,DATA,Conc,Time,n,'Obj');
        grad(i) = (or - obj)./eps_grad;
        for j = 1:length(Params)
            ej = vec;
            ej(j) = 1;
            o1 = BD_Obj_full(Params + eps_hess*ei,DATA,Conc,Time,n,'Obj');
            o2 = BD_Obj_full(Params + eps_hess*ej,DATA,Conc,Time,n,'Obj');
            o3 = BD_Obj_full(Params + eps_hess*ei + eps_hess*ej,DATA,Conc,Time,n,'Obj');
            Hess(i,j) = (o3 - o2 - o1 + obj)./(eps_hess^2);
        end
    end
end