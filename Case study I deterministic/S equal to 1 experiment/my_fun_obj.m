%% My function that return all the function value and gradient and Hessian

function f = my_fun_obj(Conc,Time,NR,theta,Observations_ave,MEAN_initial,num_sub)
    Pop = theta(1:num_sub);
    param = reshape(theta(num_sub+1:end),4,num_sub);
    f = fun(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub);
end