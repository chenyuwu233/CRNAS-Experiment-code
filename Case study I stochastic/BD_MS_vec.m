% This is a function that calculate the expected rate and variance of population at time
% t and dose concentration d given parameter Params
% Params = (r, d, b, E^n, n)
%
% Input:
%   - Params: (S x 5 matrix) that records 5 parameters of S subpopulations
%   - t: (1 x NT vector) that records NT time points
%   - d: (1 x ND vector) that records ND dosage levels
%   - n: (S x 1 vector) that records initial population of S subpopulations
%
% Output:
%   - mean: (NT x ND matrix) that records mean corresponding to each time
%   points and dosage levels
%   - var:  (NT x ND matrix) that records variance corresponding to each
%   time points and dosage levels

function [mean, var] = BD_MS_vec(Params, t,d,n)
    ND        = length(d);
    NT        = length(t);
    birth_vec = Params(:,1);
    D_effect  = log(Params(:,3)+(1-Params(:,3))./(1+(d.^Params(:,5)./Params(:,4))));
    death_mat = Params(:,2)-D_effect;
    lam    = birth_vec - death_mat;
    bd     = birth_vec + death_mat;
    mean   = zeros(NT,ND);
    var    = zeros(NT,ND);
    for i = 1:length(n)
        mean = mean + n(i)*exp(t'*lam(i,:));
        var  = var + n(i)*((bd(i,:)./lam(i,:)).*(exp(2.*t'*lam(i,:)) - exp(t'*lam(i,:))));
    end
end