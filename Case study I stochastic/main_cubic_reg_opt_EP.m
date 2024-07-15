%% Initialization
 
%% Setting the seed

parpool('local',100)
warning('off','MATLAB:integral:NonFiniteValue')
rand_seed = 13;
rng(rand_seed)

%% Initialization
    
num_optim= 20;
num_sub = 2;
lower_mixture_limit = 0.10;
initial = 1000;
max_iter = 500;
num_idx = 100;
opt_tol = 1e-6;
step_tol = 1e-6;

t_cubic_hist = [];
t_cubic_f_hist = [];
t_sqp_hist = [];
t_sqp_hess_hist = [];
t_IP_hist = [];
sol_t_cubic_hist = [];
sol_t_cubic_f_hist = [];
sol_t_sqp_hist = [];
sol_t_sqp_hess_hist = [];
sol_t_IP_hist = [];
est_cubic_hist = [];
est_cubic_f_hist = [];
est_sqp_hist = [];
est_sqp_hess_hist = [];
est_IP_hist = [];
obj_cubic_hist = [];
obj_cubic_f_hist = [];
obj_sqp_hist = [];
obj_sqp_hess_hist = [];
obj_IP_hist = [];
iter_cubic_hist = [];
iter_cubic_f_hist = [];
iter_sqp_hist = [];
iter_sqp_hess_hist = [];
iter_IP_hist = [];
modified_obj_cubic_hist = [];
modified_obj_cubic_f_hist = [];
modified_obj_sqp_hist = [];
modified_obj_sqp_hess_hist = [];
modified_obj_IP_hist = [];
good_sol_cubic_hist = [];
good_sol_cubic_f_hist = [];
good_sol_sqp_hist = [];
good_sol_sqp_hess_hist = [];
good_sol_IP_hist = [];
theta_hist = [];
true_obj_hist = [];
error_cubic_hist = [];
error_cubic_f_hist = [];
error_sqp_hist = [];
error_sqp_hess_hist = [];
error_IP_hist = [];
wrong_idx = [];


%% BD initial points
lb_BD=zeros(12,1);
ub_BD=lb_BD;
ub_BD(1)=1;ub_BD(7)=1; 
ub_BD(2)=1;ub_BD(8)=1; % natural birth rate
lb_BD(2)=0.1;lb_BD(8)=0.1;
ub_BD(3)=1;ub_BD(9)=1; % natural death rate
ub_BD(4)=1;ub_BD(10)=1; % b
lb_BD(4)=0.8;lb_BD(10)=0.8;
ub_BD(5)=10;ub_BD(11)=10; % E
lb_BD(5)=0;lb_BD(11)=0;
ub_BD(6)=10;ub_BD(12)=10; % m
lb_BD(6)=0.2;lb_BD(12)=0.2;
lb_BD   = [lb_BD;0];
ub_BD   = [ub_BD;10];
A_BD    = [0,1,-1,0,0,0,0,0,0,0,0,0,0;
           0,0,0,0,0,0,0,1,-1,0,0,0,0;
           0,-1,1,0,0,0,0,0,0,0,0,0,0;
           0,0,0,0,0,0,0,-1,1,0,0,0,0];
b_BD    = [0.1;0.1;0;0];
A_eq    = [1,0,0,0,0,0,1,0,0,0,0,0,0];
b_eq    = 1;


%% DATA generation
lb_GE = lb_BD;
ub_GE = ub_BD;
lb_GE(5)  = 0.05;
ub_GE(5)  = 0.1;
lb_GE(11) = 0.5;
ub_GE(11) = 2.5;
lb_GE(1)  = 0.3;
ub_GE(1)  = 0.5;
lb_GE(4)  = 0.8;
ub_GE(4)  = 1;
lb_GE(10)  = 0.8;
ub_GE(10)  = 1;
ub_GE(6)  = 5;
lb_GE(6)  = 1.5;
ub_GE(12) = 5;
lb_GE(12) = 1.5;

NT      = 13; 
NC      = 11;
NR      = 13;



Conc = 10^(6)*[0 31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6)];
Time =   3*(0:NT-1);

%% Change the lb and ub

    % lb_BD_opt = lb_BD;
    % ub_BD_opt = ub_BD;

    ub_BD_sqp = [1,1,1,1,Inf,Inf,1,1,1,1,Inf,Inf,10];
    ub_BD_cubic = [1,1,1,1,Inf,Inf,1,1,1,1,Inf,Inf,10];
    lb_BD_sqp = [0,0,0,0,0,0.01,0,0,0,0,0,0.01,0];
    lb_BD_cubic = [0,0,0,0,0,0.01,0,0,0,0,0,0.01,0];

   

    %%
parfor e_num = 1:num_idx

    
    
        x_init_BD    = [];
        for n=1:num_optim
            x02      = rand(length(ub_BD),1).*(ub_BD-lb_BD) + lb_BD;
            x02(3) = max(0,x02(2) - rand*0.1);
            x02(9) = max(x02(8)  - rand*0.1,0);
            x02(7) = 1-x02(1); 
            x_init_BD   = [x_init_BD,x02];
        end




%%

        theta = rand(13,1).*(ub_GE-lb_GE) + lb_GE;
        theta(3) = max(0,theta(2) - rand * 0.1);
        theta(9) = max(theta(8) - rand * 0.1,0);
        theta(7) = 1 - theta(1);
        theta    = GR_sort(theta,max(Conc));
        
        Prec_p   = [];
        Prec_GR1 = [];
        Prec_GR2 = [];
        indi_ip  = get_indi(theta,Conc(end));
%%
        DATA = sto_gen_bd_EP(NR,Conc,Time,initial,theta);
    
        


    


    %% Change of variable

    x_init = x_init_BD;
    theta_c = theta;
    % x_init(5,:) = x_init(5,:).^x_init(6,:);
    % x_init(11,:) = x_init(11,:).^x_init(12,:);
    theta_c(5,:) = theta_c(5,:).^theta_c(6,:);
    theta_c(11,:) = theta_c(11,:).^theta_c(12,:);

    %% Check for objective

    
    DATA_old = permute(DATA,[3,2,1]);
    
    tic
    obj_old = BD_Obj_full(theta_c,DATA,Conc,Time,initial,'Obj');
    t1 = toc
    tic
    [obj_new,grad,Hess]= BD_Obj_full(theta_c,DATA,Conc,Time,initial,'Hess');
    t2 = toc
    tic
    [obj_new2,grad2]= BD_Obj_grad(theta_c,DATA,Conc,Time,initial,'grad');
    t2_grad = toc
    tic
    [obj_new_f,grad_f,Hess_f]= BD_Obj_full_finite(theta_c,DATA,Conc,Time,initial,1e-10,1e-7);
    t3 = toc

    
% pause

    
    %% Cubic reg

    val_cubic=[];
    nor_dis_i=[];   
    nor_dis_i_2=[];
    x_list_cubic=[];
    idx=[];
    iter_cubic = [];
    sol_t_cubic = [];



    % tic
        for ii=1:num_optim
            tic
            [param,val,nor_dis,nor_dis2,iter]=mysolve_EP(Conc,Time,NR,x_init(:,ii),DATA,initial,2,lb_BD_cubic',ub_BD_cubic',max_iter,opt_tol,step_tol);
            % if val(end) <= obj_old
                val_cubic=[val_cubic,val(end)];
                param(5) = param(5)^(1/param(6));
                param(11) = param(11)^(1/param(12));
                param = GR_sort(param,max(Conc));
                iter_cubic = [iter_cubic,iter];
                x_list_cubic=[x_list_cubic,param];
                idx=[idx,ii];
                tii = toc;
                sol_t_cubic = [sol_t_cubic,tii];
            % end
        end

    % t_cubic=toc

        %% Cubic reg finite

    val_cubic_f=[];
    nor_dis_i=[];   
    nor_dis_i_2=[];
    x_list_cubic_f=[];
    idx=[];
    iter_cubic_f = [];
    sol_t_cubic_f = [];
    



    % tic
        for ii=1:num_optim
            tic
            [param,val,nor_dis,nor_dis2,iter]=mysolve_EP_f(Conc,Time,NR,x_init(:,ii),DATA,initial,2,lb_BD_cubic',ub_BD_cubic',max_iter,opt_tol,step_tol);
            % if val(end)<=obj_old
                val_cubic_f=[val_cubic_f,val(end)];
                param(5) = param(5)^(1/param(6));
                param(11) = param(11)^(1/param(12));
                param = GR_sort(param,max(Conc));
                iter_cubic_f = [iter_cubic_f,iter];
                x_list_cubic_f=[x_list_cubic_f,param];
                idx=[idx,ii];
                tii = toc;
                sol_t_cubic_f = [sol_t_cubic_f,tii];
            % end
        end

    % t_cubic_f=toc


    %%  Fmincon SQP without Hessian


        options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
            ,'algorithm','sqp','OptimalityTolerance',opt_tol,'StepTolerance',step_tol);






    % tic


    val_sqp  = [];
    params2_sqp = [];
    grad2_sqp = [];
    iter_sqp = [];
    sol_t_sqp = [];


%     func_sqp = @(x) BD_Obj(x,DATA_old,Conc,Time,NR,2);
    func_sqp = @(x) BD_Obj_full(x,DATA,Conc,Time,initial,'Obj');

    for ii=1:num_optim
        tic
        [xx_sqp,ff_sqp,ef_sqp,out_sqp,~,g_sqp,~]  = fmincon(func_sqp,x_init(:,ii),[],[],A_eq,b_eq,lb_BD_sqp,ub_BD_sqp,[],options1);
        % if ff_sqp <= obj_old
            xx_sqp(5) = xx_sqp(5)^(1/xx_sqp(6));
            xx_sqp(11) = xx_sqp(11)^(1/xx_sqp(12));
            xx_sqp    = GR_sort(xx_sqp,max(Conc));
            xx_sqp = GR_sort(xx_sqp,max(Conc));
            val_sqp    = [val_sqp, ff_sqp];
            params2_sqp  = [params2_sqp,xx_sqp];
            iter_sqp = [iter_sqp,out_sqp.iterations];
            tii = toc;
            sol_t_sqp = [sol_t_sqp,tii];
            % grad2_sqp    = [grad2_sqp,g_sqp];
        % end
    end
    % t_sqp=toc



    %%  Fmincon SQP with Hessian



    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
        ,'algorithm','sqp','SpecifyObjectiveGradient',true,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol);

    % tic


    val_sqp_hess  = [];
    params2_sqp_hess = [];
    iter_sqp_hess = [];
    sol_t_sqp_hess = [];
    func_sqp = @(x) BD_Obj_grad(x,DATA,Conc,Time,initial,'grad');

    for ii=1:num_optim
        tic
        [xx_sqp,ff_sqp,ef_sqp,out_sqp,~,g_sqp,~]  = fmincon(func_sqp,x_init(:,ii),[],[],A_eq,b_eq,lb_BD_sqp,ub_BD_sqp,[],options1);
        % if ff_sqp <= obj_old
            xx_sqp(5) = xx_sqp(5)^(1/xx_sqp(6));
            xx_sqp(11) = xx_sqp(11)^(1/xx_sqp(12));
            xx_sqp    = GR_sort(xx_sqp,max(Conc));
    
            val_sqp_hess    = [val_sqp_hess, ff_sqp];
            params2_sqp_hess  = [params2_sqp_hess,xx_sqp];
            iter_sqp_hess     = [iter_sqp_hess,out_sqp.iterations];
            tii = toc;
            sol_t_sqp_hess    = [sol_t_sqp_hess,tii]; 
        % end
    end
    % t_sqp_hess=toc


%%  Fmincon IP with Hessian



    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
        ,'SpecifyObjectiveGradient',true,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol);

    % tic


    val_IP  = [];
    params2_IP = [];
    func_IP = @(x) BD_Obj_full(x,DATA,Conc,Time,initial,'Hess');
    g_IP_hist = [];
    ef_IP_hist = [];
    iter_IP = [];
    sol_t_IP = [];
    for ii=1:num_optim
        tic
        [xx_IP,ff_IP,ef_IP,out_IP,~,g_IP,~]  = fmincon(func_IP,x_init(:,ii),[],[],A_eq,b_eq,lb_BD_sqp,ub_BD_sqp,[],options1);
        % if ff_IP <= obj_old
            xx_IP(5) = xx_IP(5)^(1/xx_IP(6));
            xx_IP(11) = xx_IP(11)^(1/xx_IP(12));
            xx_IP    = GR_sort(xx_IP,max(Conc));
    
            val_IP    = [val_IP, ff_IP];
            params2_IP  = [params2_IP,xx_IP];
            ef_IP_hist  = [ef_IP_hist,ef_IP];
            g_IP_hist   = [g_IP_hist,g_IP];
            iter_IP     = [iter_IP,out_IP.iterations];
            tii = toc;
            sol_t_IP    = [sol_t_IP,tii];
        % end
    end
    % t_IP=toc







    %% Recording
    try
        num_best=1;
        [B,I]=mink(val_cubic,num_best);
        [B_f,I_f] = mink(val_cubic_f,num_best);
        [B_sqp,I_sqp]=mink(val_sqp,num_best);
        [B_sqp_hess,I_sqp_hess]=mink(val_sqp_hess,num_best);
        [B_IP,I_IP]=mink(val_IP,num_best);
        % [~,I]=min(iter_cubic);
        % [~,I_f]=min(iter_cubic_f);
        % [~,I_sqp]=min(iter_sqp);
        % [~,I_sqp_hess]=min(iter_sqp_hess);
        % [~,I_IP]=min(iter_IP);
        % B = val_cubic(I);
        % B_f = val_cubic_f(I_f);
        % B_sqp = val_sqp(I_sqp);
        % B_sqp_hess = val_sqp_hess(I_sqp_hess);
        % B_IP = val_IP(I_IP);


    
        theta_hist = [theta_hist,theta];
        true_obj_hist = [true_obj_hist,obj_new];
        t_cubic_hist = [t_cubic_hist,sum(sol_t_cubic)];
        t_cubic_f_hist = [t_cubic_f_hist,sum(sol_t_cubic_f)];
        t_sqp_hist = [t_sqp_hist,sum(sol_t_sqp)];
        t_sqp_hess_hist = [t_sqp_hess_hist,sum(sol_t_sqp_hess)];
        t_IP_hist = [t_IP_hist,sum(sol_t_IP)];
        est_cubic_hist = [est_cubic_hist,x_list_cubic(:,I)];
        est_cubic_f_hist = [est_cubic_f_hist,x_list_cubic_f(:,I_f)];
        est_sqp_hist = [est_sqp_hist,params2_sqp(:,I_sqp)];
        est_sqp_hess_hist = [est_sqp_hess_hist, params2_sqp_hess(:,I_sqp_hess)];
        est_IP_hist = [est_IP_hist,params2_IP(:,I_IP)];
        obj_cubic_hist = [obj_cubic_hist,B];
        obj_cubic_f_hist = [obj_cubic_f_hist,B_f];
        obj_sqp_hist = [obj_sqp_hist,B_sqp];
        obj_sqp_hess_hist = [obj_sqp_hess_hist,B_sqp_hess];
        obj_IP_hist = [obj_IP_hist,B_IP];
        error_cubic_hist = [error_cubic_hist,abs(x_list_cubic(:,I) - theta)];
        error_cubic_f_hist = [error_cubic_f_hist,abs(x_list_cubic_f(:,I_f) - theta)];
        error_sqp_hist = [error_sqp_hist,abs(params2_sqp(:,I_sqp) - theta)];
        error_sqp_hess_hist = [error_sqp_hess_hist,abs(params2_sqp_hess(:,I_sqp_hess) - theta)];
        error_IP_hist = [error_IP_hist,abs(params2_IP(:,I_IP) - theta)];
        iter_cubic_hist = [iter_cubic_hist,iter_cubic(I)];
        iter_cubic_f_hist = [iter_cubic_f_hist,iter_cubic_f(I_f)];
        iter_sqp_hist = [iter_sqp_hist,iter_sqp(I_sqp)];
        iter_sqp_hess_hist = [iter_sqp_hess_hist,iter_sqp_hess(I_sqp_hess)];
        iter_IP_hist = [iter_IP_hist,iter_IP(I_IP)];
        modified_obj_cubic_hist = [modified_obj_cubic_hist,B./obj_old];
        modified_obj_cubic_f_hist = [modified_obj_cubic_f_hist,B_f./obj_old];
        modified_obj_sqp_hist = [modified_obj_sqp_hist,B_sqp./obj_old];
        modified_obj_sqp_hess_hist = [modified_obj_sqp_hess_hist,B_sqp_hess./obj_old];
        modified_obj_IP_hist = [modified_obj_IP_hist,B_IP./obj_old];
        sol_t_cubic_hist = [sol_t_cubic_hist,sol_t_cubic(I)];
        sol_t_cubic_f_hist = [sol_t_cubic_f_hist,sol_t_cubic_f(I_f)];
        sol_t_sqp_hist = [sol_t_sqp_hist,sol_t_sqp(I_sqp)];
        sol_t_sqp_hess_hist = [sol_t_sqp_hess_hist,sol_t_sqp_hess(I_sqp_hess)];
        sol_t_IP_hist = [sol_t_IP_hist,sol_t_IP(I_IP)];
        good_sol_cubic_hist = [good_sol_cubic_hist,sum(val_cubic<=obj_old)];
        good_sol_cubic_f_hist = [good_sol_cubic_f_hist,sum(val_cubic_f<=obj_old)];
        good_sol_sqp_hist = [good_sol_sqp_hist,sum(val_sqp<=obj_old)];
        good_sol_sqp_hess_hist = [good_sol_sqp_hess_hist,sum(val_sqp_hess<=obj_old)];
        good_sol_IP_hist = [good_sol_IP_hist,sum(val_IP<=obj_old)];
    catch
        wi = 0
        if isempty(val_cubic) 
            wi = wi+1;
        elseif isempty(val_cubic_f) 
            wi = wi+2;
        elseif isempty(val_sqp) 
            wi = wi+4;
        elseif isempty(val_sqp_hess) 
            wi = wi+8;
        elseif isempty(val_IP)
            wi = wi+16;
        end
        wrong_idx = [wrong_idx,wi];
    end
end

name = strcat('Result/BD_opt',num2str(num_optim),'_',num2str(num_idx),'_uc.mat');
save(name)

