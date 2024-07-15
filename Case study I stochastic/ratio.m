% This is a function that calculate the gradient of ratio function such
% that ratio(beta,nu,H) = (beta+nu-H)/(beta-nu+H)
%
% Input:
%   - d: (scalar) dosage levels
%   - beta: (scalar)
%   - nu: (scalar)
%   - b: (scalar) 
%   - E: (scalar) changed E = old E ^ m;
%   - m: (scalar) 
%
% Output:
%   - obj: (scalar)
%   - grad: (NT x ND x 5) [p,beta,nu,b,E,m]
%   - Hess: (5 x 5 matrix) [p,beta,nu,b,E,m]
%
% Requirment:
%   - Hill(d,b,E,m)

function [obj,grad,Hess] = ratio(d,beta,nu,b,E,m)
    [obj_H,grad_H,Hess_H] = Hill(d,b,E,m);
    bn  = (beta + nu - obj_H);
    lam = (beta - nu + obj_H);
    obj = bn/lam;
    grad = [2*(obj_H - nu)/lam^2, 2*beta/lam^2, -lam^(-1).*grad_H - bn/lam^2.*grad_H];
    Hess = zeros(5,5);
    Hess(3:5,3:5) = lam^(-2).* grad_H'*grad_H - lam^(-1).* Hess_H + (lam+2*bn)/lam^3.*grad_H'*grad_H - bn/lam^2.* Hess_H;
    Hess(1,1) = -4*(obj_H-nu)/lam^(3);
    Hess(2,2) = 4*beta/lam^(3);
    Hess(1,2) = (-2*lam + 4*(obj_H-nu))./lam^(3);
    Hess(2,1) = Hess(1,2);
    Hess(1,3:5) = (2*lam-4*(obj_H-nu))/lam^(3).*grad_H;
    Hess(3:5,1) = (2*lam-4*(obj_H-nu))/lam^(3).*grad_H;
    Hess(3:5,2) = -4*beta/lam^(3).*grad_H;
    Hess(2,3:5) = -4*beta/lam^(3).*grad_H;
end