%% My function that return all the function value and gradient and Hessian

function [f,g] = my_fun_mono_grad(Conc,Time,NR,theta,Observations_ave,MEAN_initial)
    % Pop = theta(1:num_sub);
    param = theta;
    f = fun(Conc,Time,NR,1,param,Observations_ave,MEAN_initial,1);
    [~,g_param] = grad(Conc,Time,NR,1,param,Observations_ave,MEAN_initial,1);
    g = [vec(g_param)];
end