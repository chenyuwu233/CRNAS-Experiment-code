%% Function to reshape the param where param: [Initial_pop, param_sub, noise]

function [Pop_vec, Param_mat, noise_vec] = reshape_param(param,num_sub)

    noise_vec = param(end-1:end);
    param(end-1:end) = [];
    Pop_vec = param(1:num_sub-1);
    Pop_vec = [Pop_vec, 1-sum(Pop_vec)];
    Param_mat = reshape(param(num_sub:end),4,num_sub);
    
end