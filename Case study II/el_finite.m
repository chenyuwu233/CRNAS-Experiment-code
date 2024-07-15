%%
% params = [a,b,n,rho,g,m,d,s]


function [obj,grad,Hess] = el_finite(Params,DATA,Time,Init,eps_grad,eps_hess,varargin)

    obj = el(Params,DATA,Time,Init);
    grad = zeros(length(Params),1);
    Hess = zeros(length(Params),length(Params));
    if nargin == 7
        vec  = zeros(size(Params));
        for i = 1:length(Params)
            ei = vec;
            ei(i) = 1;
            ol = el(Params - eps_grad*ei,DATA,Time,Init);
            or = el(Params + eps_grad*ei,DATA,Time,Init);
            grad(i) = (or - ol)./(2*eps_grad);
            for j = 1:length(Params)
                ej = vec;
                ej(j) = 1;
                o1 = el(Params + eps_hess*ei,DATA,Time,Init);
                o2 = el(Params + eps_hess*ej,DATA,Time,Init);
                o3 = el(Params + eps_hess*ei + eps_hess*ej,DATA,Time,Init);
                Hess(i,j) = (o3 - o2 - o1 + obj)./(eps_hess^2);
            end
        end
    end
end