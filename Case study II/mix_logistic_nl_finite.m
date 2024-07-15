function [obj,grad,Hess] = mix_logistic_nl_finite(theta,DATA,cv,init,eps_grad,eps_hess,varargin)


    for i = 1:length(DATA)

        func_i = @(x) mix_logistic_nl_obj(x,cv(i),DATA(i),init);
        obj = func_i(theta);
        grad = zeros(length(theta),1);
        Hess = zeros(length(theta),length(theta));
    
        if nargin == 7
            vec  = zeros(size(theta));
            for j = 1:length(theta)
                ei = vec;
                ei(j) = 1;
                ol = func_i(theta - eps_grad*ei);
                or = func_i(theta + eps_grad*ei);
                grad(j) = (or - ol)./(2*eps_grad);
                for k = 1:length(theta)
                    ej = vec;
                    ej(k) = 1;
                    o1 = func_i(theta + eps_hess*ei);
                    o2 = func_i(theta + eps_hess*ej);
                    o3 = func_i(theta + eps_hess*ei + eps_hess*ej);
                    Hess(j,k) = (o3 - o2 - o1 + obj)./(eps_hess^2);
                end
            end
        end

    end
end