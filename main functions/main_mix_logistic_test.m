addpath('mix_logistic_test')

rng(13)
%%
ND = 20;
% num_sub = 2;
num_opt = 20;
num_idx = 100;
max_iter = 500;
opt_tol = 1e-6;
step_tol = 1e-6;

%% Test for logistic growth estimation
gen_lb_1 = [0,0,0];
gen_ub_1 = [1,1,1];

gen_lb_2 = [0,2,2];
gen_ub_2 = [1,3,3];

theta_1 = rand(1,3).*(gen_ub_1-gen_lb_1)+gen_lb_1;
theta_2 = rand(1,3).*(gen_ub_2-gen_lb_2)+gen_lb_2;
theta_2(1) = 1 - theta_1(1);
theta = [theta_1,theta_2];
init = 1000;


Dx = linspace(0,10,ND);

DATA = [];
for i = 1:ND
    est  = init*theta_1(1)*logistic(Dx(i),1,theta_1(2),theta_1(3)) ... 
          +init*theta_2(1)*logistic(Dx(i),1,theta_2(2),theta_2(3));
    DATA = [DATA,est];
end
    
plot(Dx,DATA)

%% Initialize points
    
lb  = [0,0,0,0,0,0];
ub  = [1,10,10,1,10,10];
Aeq = [1,0,0,1,0,0];
beq = 1;

num_opt = 20;
x_init = [];
for i = 1:num_opt
    xi = rand(1,6).*ub;
    xi(4) = 1-xi(1);
    x_init = [x_init;xi];     
end

%% Cubic initialization

func_full = @mix_logistic_nl_full;
func_obj  = @test_obj;
func_inputs.init = init;
func_inputs.num_sub = 2;
func_inputs.cmd = 2;


%% Optimization

val_cubic_f=[];
x_list_cubic_f=[];
iter_hist_cubic = [];
idx=[];
tic
for ii=1:num_opt
    [param,val,iter]=CRNAS(func_full,func_obj,func_inputs,x_init(ii,:)',DATA,Dx,Aeq,beq,lb',ub',max_iter,opt_tol,step_tol);
    val_cubic_f=[val_cubic_f,val(end)];
    x_list_cubic_f=[x_list_cubic_f,param];
    iter_hist_cubic = [iter_hist_cubic,iter];
    idx=[idx,ii];
end
t_cubic=toc
[fval_cubic,oi_cubic] = min(val_cubic_f);
opt_xx_pe_cubic = x_list_cubic_f(:,oi_cubic);
iter_cubic = iter_hist_cubic(oi_cubic);