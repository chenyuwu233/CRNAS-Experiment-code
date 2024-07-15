% This is a function that calculate the gradient of the Hill(changed) at specific
% variable, H(d;b,E,m) = log(b + (1-b)/(1+(d^m/E)));
%
% Input:
%   - d: (scalar) dosage level
%   - b: (scalar) 
%   - E: (scalar) changed E = old E ^ m;
%   - m: (scalar) 
%
% Output:
%   - obj: (scalar)
%   - grad: (1 x 3 vector) [b,E,m]
%   - Hess: (3 x 3 matrix) [b,E,m]
%
% Requirment:
%   - grad_pop_H(vec,d,cmd)

function [obj,grad,Hess] = Hill(d,b,E,m)
    %% Objective
    obj = log(b + (1-b)/(1+d^m/E));
    vec = [0,b,E,m];
    [gb,gE,gm,Hess] = grad_pop_H(vec,d,1);
    grad = [gb,gE,gm];

    %% Gradient
%     grad_share = 1/(b + (1-b)/(1+d^m/E)); % Log part
%     grad = zeros(1,3);
%     grad(1) = grad_share*(1 - 1/(1+d^m/E));
%     grad(2) = grad_share*(1-b)*(1+d^m/E)^(-2)*(d^m/E^2);
%     grad(3) = grad_share*(1-b)*(1+d^m/E)^(-2)*(-d^m*log(d)/E);

    %% Hessian
%     Hess_share = -(b + (1-b)/(1+d^m/E))^(-2);
%     Hess = zeros(3,3);
%     Hess(1,1) = Hess_share*(1 - 1/(1+d^m/E))*(1 - 1/(1+d^m/E));
%     Hess(2,2) = Hess_share*(1-b)*(1+d^m/E)^(-2)*(d^m/E^2)*(1-b)*(1+d^m/E)^(-2)*(d^m/E^2) + grad(2)*(2*(1-b)*(1+d^m/E)^(-3)*(d^m/E^2)^2 - 2*(1-b)*(1+d^m/E)^(-2)*(d^m/E^3));
%     Hess(3,3) = Hess_share*(1-b)*(1+d^m/E)^(-2)*(-d^m*log(d)/E)*(1-b)*(1+d^m/E)^(-2)*(-d^m*log(d)/E) + grad(3)*(2*(1-b)*(1+d^m/E)^(-3)*(d^m*log(d)/E)^2 - (1-b)*(1+d^m/E)^(-2)*(d^m*log(d)*log(d)/E));
%     Hess(1,2) = Hess_share*(1 - 1/(1+d^m/E))*(1-b)*(1+d^m/E)^(-2)*(d^m/E^2) + grad_share*(-(1+d^m/E)^(-2)*d^m/E^2);
%     Hess(2,1) = Hess(1,2);
%     Hess(1,3) = Hess_share*(1 - 1/(1+d^m/E))*(1-b)*(1+d^m/E)^(-2)*(-d^m*log(d)/E) + grad_share*((1+d^m/E)^(-2)*d^m*log(d)/E);
%     Hess(3,1) = Hess(1,3);
%     Hess(2,3) = Hess_share*(1-b)*(1+d^m/E)^(-2)*(d^m/E^2)*(1-b)*(1+d^m/E)^(-2)*(-d^m*log(d)/E) + grad_share*(-2*(1-b)*(1+d^m/E)^(-3)*d^(2*m)*log(d)/E^3 + (1-b)*(1+d^m/E)^(-2)*d^m*log(d)/E^2);
%     Hess(3,2) = Hess(2,3);


end