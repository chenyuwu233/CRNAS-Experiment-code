% %% Initialization
% 
% num_optim= 5;
% num_sub = 2;
% lower_mixture_limit = 0.10;
% initial = 1000;
% max_iter = 150;
% num_idx = 100;
% opt_tol = 1e-2;
% step_tol = 1e-4;
% 
% t_cubic_hist = [];
% t_cubic_f_hist = [];
% t_sqp_hist = [];
% t_sqp_hess_hist = [];
% t_IP_hist = [];
% est_cubic_hist = [];
% est_cubic_f_hist = [];
% est_sqp_hist = [];
% est_sqp_hess_hist = [];
% est_IP_hist = [];
% obj_cubic_hist = [];
% obj_cubic_f_hist = [];
% obj_sqp_hist = [];
% obj_sqp_hess_hist = [];
% obj_IP_hist = [];
% theta_hist = [];
% true_obj_hist = [];
% error_cubic_hist = [];
% error_cubic_f_hist = [];
% error_sqp_hist = [];
% error_sqp_hess_hist = [];
% error_IP_hist = [];
% 
% 
% %% BD initial points
% lb_BD=zeros(12,1);
% ub_BD=lb_BD;
% ub_BD(1)=1;ub_BD(7)=1; 
% ub_BD(2)=1;ub_BD(8)=1; % natural birth rate
% lb_BD(2)=0.1;lb_BD(8)=0.1;
% ub_BD(3)=1;ub_BD(9)=1; % natural death rate
% ub_BD(4)=1;ub_BD(10)=1; % b
% lb_BD(4)=0.8;lb_BD(10)=0.8;
% ub_BD(5)=10;ub_BD(11)=10; % E
% lb_BD(5)=0;lb_BD(11)=0;
% ub_BD(6)=10;ub_BD(12)=10; % m
% lb_BD(6)=0.2;lb_BD(12)=0.2;
% lb_BD   = [lb_BD;0];
% ub_BD   = [ub_BD;10];
% A_BD    = [0,1,-1,0,0,0,0,0,0,0,0,0,0;
%            0,0,0,0,0,0,0,1,-1,0,0,0,0;
%            0,-1,1,0,0,0,0,0,0,0,0,0,0;
%            0,0,0,0,0,0,0,-1,1,0,0,0,0];
% b_BD    = [0.1;0.1;0;0];
% A_eq    = [1,0,0,0,0,0,1,0,0,0,0,0,0];
% b_eq    = 1;
% 
% 
% %% DATA generation
% lb_GE = lb_BD;
% ub_GE = ub_BD;
% lb_GE(5)  = 0.05;
% ub_GE(5)  = 0.1;
% lb_GE(11) = 0.5;
% ub_GE(11) = 2.5;
% lb_GE(1)  = 0;
% lb_GE(4)  = 0.8;
% ub_GE(4)  = 1;
% lb_GE(10)  = 0.8;
% ub_GE(10)  = 1;
% ub_GE(6)  = 5;
% lb_GE(6)  = 1.5;
% ub_GE(12) = 5;
% lb_GE(12) = 1.5;
% 
% NT      = 13; 
% NC      = 11;
% NR      = 13;
% 
% 
% 
% Conc = 10^(6)*[0 31.25*10^(-9) 62.5*10^(-9) 125*10^(-9) 250*10^(-9) 375*10^(-9) 500*10^(-9) 1.25*10^(-6) 2.5*10^(-6) 3.75*10^(-6) 5*10^(-6)];
% Time =   3*(0:NT-1);
% 
% %% Change the lb and ub
% 
% % lb_BD_opt = lb_BD;
% % ub_BD_opt = ub_BD;
% 
% ub_BD_sqp = [1,1,1,1,Inf,Inf,1,1,1,1,Inf,Inf,10];
% ub_BD_cubic = [1,1,1,1,Inf,Inf,1,1,1,1,Inf,Inf,10];
% lb_BD_sqp = [0,0,0,0,0,0.01,0,0,0,0,0,0.01,0];
% lb_BD_cubic = [0,0,0,0,0,0.01,0,0,0,0,0,0.01,0];
% 
% 
% 
% 
% 
%     x_init_BD    = [];
%     for n=1:num_optim
%         x02      = rand(length(ub_BD),1).*(ub_BD-lb_BD) + lb_BD;
%         x02(3) = max(0,x02(2) - rand*0.1);
%         x02(9) = max(x02(8)  - rand*0.1,0);
%         x02(7) = 1-x02(1); 
%         x_init_BD   = [x_init_BD,x02];
%     end
% 
% 
% 
% 
% %%
% 
%     theta = rand(13,1).*(ub_GE-lb_GE) + lb_GE;
%     theta(3) = max(0,theta(2) - rand * 0.1);
%     theta(9) = max(theta(8) - rand * 0.1,0);
%     theta(7) = 1 - theta(1);
%     theta    = GR_sort(theta,max(Conc));
% 
%     Prec_p   = [];
%     Prec_GR1 = [];
%     Prec_GR2 = [];
%     indi_ip  = get_indi(theta,Conc(end));
% 
%     DATA = sto_gen_bd_EP(NR,Conc,Time,initial,theta);
% 
% 
% 
% 
% 
% 
% 
% %% Change of variable
% 
% x_init = x_init_BD;
% theta_c = theta;
% % x_init(5,:) = x_init(5,:).^x_init(6,:);
% % x_init(11,:) = x_init(11,:).^x_init(12,:);
% theta_c(5,:) = theta_c(5,:).^theta_c(6,:);
% theta_c(11,:) = theta_c(11,:).^theta_c(12,:);



    % Cubic reg

val_cubic=[];
nor_dis_i=[];   
nor_dis_i_2=[];
x_list_cubic=[];
idx=[];




tic
    for ii=1:num_optim
        [param,val,nor_dis,nor_dis2]=mysolve_EP(Conc,Time,NR,x_init(:,ii),DATA,initial,2,lb_BD_cubic',ub_BD_cubic',max_iter,opt_tol,step_tol);
        val_cubic=[val_cubic,val(end)];
        param(5) = param(5)^(1/param(6));
        param(11) = param(11)^(1/param(12));
        param = GR_sort(param,max(Conc));
        x_list_cubic=[x_list_cubic,param];
        idx=[idx,ii];
    end

t_cubic=toc



% %%  Fmincon IP with Hessian
% 
% 
% 
%     options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
%         ,'SpecifyObjectiveGradient',true,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol);
% 
%     tic
% 
% 
%     fval2_IP  = [];
%     params2_IP = [];
%     func_IP = @(x) BD_Obj_full(x,DATA,Conc,Time,initial,'Hess');
%     g_IP_hist = [];
%     ef_IP_hist = [];
%     for ii=1:num_optim
%         [xx_IP,ff_IP,ef_IP,out_IP,~,g_IP,~]  = fmincon(func_IP,x_init(:,ii),A_BD,b_BD,A_eq,b_eq,lb_BD_sqp,ub_BD_sqp,[],options1);
%         xx_IP(5) = xx_IP(5)^(1/xx_IP(6));
%         xx_IP(11) = xx_IP(11)^(1/xx_IP(12));
%         xx_IP    = GR_sort(xx_IP,max(Conc));
% 
%         fval2_IP    = [fval2_IP, ff_IP];
%         params2_IP  = [params2_IP,xx_IP];
%         ef_IP_hist  = [ef_IP_hist,ef_IP];
%         g_IP_hist   = [g_IP_hist,g_IP];
%     end
%     t_IP=toc
