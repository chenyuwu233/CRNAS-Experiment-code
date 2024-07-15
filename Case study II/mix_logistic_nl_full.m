%% General skeleton 

function [obj,grad,Hess] = mix_logistic_nl_full(theta,DATA,cv,init,num_sub,cmd)
    if cmd >= 0
        obj = 0;
        for i = 1:length(DATA)
            obj  = obj  + mix_logistic_nl_obj(theta,cv(i),DATA(i),init);
        end
    end
    if cmd>=1
        grad = zeros(3*num_sub,1);
        for i = 1:length(DATA)
            grad = grad + mix_logistic_nl_grad(theta,cv(i),DATA(i),init);
        end
    end
    if cmd>=2
        Hess = zeros(3*num_sub,3*num_sub);
        for i = 1:length(DATA)
            Hess = Hess + mix_logistic_nl_Hess(theta,cv(i),DATA(i),init);
        end
    end    
end