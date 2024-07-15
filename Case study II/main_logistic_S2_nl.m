rng(3)

%% Quantitative experiment

cubic_val_hist = [];
cubic_err_hist = [];
cubic_t_hist   = [];
cubic_est_hist = [];
cubic_iter_hist = [];
cubic_number_guess_hist = [];
fmin_IP_val_hist  = [];
fmin_IP_err_hist  = [];
fmin_IP_t_hist    = [];
fmin_IP_est_hist  = [];
fmin_IP_iter_hist = [];
fmin_IP_number_guess_hist = [];
fmin_IP_hess_val_hist  = [];
fmin_IP_hess_err_hist  = [];
fmin_IP_hess_t_hist    = [];
fmin_IP_hess_est_hist  = [];
fmin_IP_hess_iter_hist = [];
fmin_IP_hess_number_guess_hist = [];
fmin_sqp_val_hist  = [];
fmin_sqp_err_hist  = [];
fmin_sqp_t_hist    = [];
fmin_sqp_est_hist  = [];
fmin_sqp_iter_hist = [];
fmin_sqp_number_guess_hist = [];
fmin_sqp_grad_val_hist  = [];
fmin_sqp_grad_err_hist  = [];
fmin_sqp_grad_t_hist    = [];
fmin_sqp_grad_est_hist  = [];
fmin_sqp_grad_iter_hist = [];
fmin_sqp_grad_number_guess_hist = [];
theta_hist = [];
ND = 20;
% num_sub = 2;
num_opt = 20;
num_idx = 100;
max_iter = 500;
opt_tol = 1e-6;
step_tol = 1e-6;




%%

for ii = 1:num_idx


%% Test for logistic growth estimation


    gen_lb_1 = [0,0,0];
    gen_ub_1 = [1,1,1];

    gen_lb_2 = [0,2,2];
    gen_ub_2 = [1,3,3];

    theta_1 = rand(1,3).*(gen_ub_1-gen_lb_1)+gen_lb_1;
    theta_2 = rand(1,3).*(gen_ub_2-gen_lb_2)+gen_lb_2;
    % theta_1(3) = theta_1(2) * theta_1(3);
    % theta_2(3) = theta_2(2) * theta_2(3);
    theta_2(1) = 1 - theta_1(1);
    theta = [theta_1,theta_2];
    init = 1000;

    theta_hist = [theta_hist;theta];
    
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
    

%%
    ttheta = theta;
    % mix_obj = @(x) mix_logistic_full(x,DATA,Dx,init,2,0);

    [obj,grad,Hess] = mix_logistic_nl_full(ttheta,DATA,Dx,init,2,2)
    
    eig(Hess)

    [obj_f,grad_f,Hess_f] = mix_logistic_nl_finite(ttheta,DATA,Dx,init,1e-9,1e-9,1)


    eig(Hess_f) 

    norm(Hess - Hess_f)
    % eps = 1e-5;
    % 
    % g_eps = zeros(8,1);
    % H_eps = zeros(8,8);
    % 
    % for i = 1:8
    %     eps_vec = zeros(1,8);
    %     eps_vec(i) = eps;
    %     ttheta_l = ttheta - eps_vec;
    %     ttheta_r = ttheta + eps_vec;
    %     g_eps(i) = (mix_obj(ttheta_r) - mix_obj(ttheta_l))./(2*eps);
    % end
    % 
    % for i = 1:8
    %     tmp_i = zeros(1,8);
    %     tmp_i(i) = eps;
    %     for j = 1:8
    %         tmp_j = zeros(1,8);
    %         tmp_j(j) = eps;
    %         H_eps(i,j) = (mix_obj(ttheta+tmp_i+tmp_j) - mix_obj(ttheta+tmp_i)-mix_obj(ttheta+tmp_j)+mix_obj(ttheta))./(eps^2);
    %     end
    % end
    % 
    % norm(H_eps-Hess)
    % 
    % 
    % 


    
    %% fmincon IP optimize
    func = @(x) mix_logistic_nl_full(x,DATA,Dx,init,2,0);
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
            ,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol);
    
    fval_hist_pe_IP   = [];
    params_hist_pe_IP = [];
    iter_hist_IP      = [];
    tic
    for j = 1:num_opt
        [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],Aeq,beq,lb,ub,[],options1); 
        fval_hist_pe_IP = [fval_hist_pe_IP,ff];
        params_hist_pe_IP = [params_hist_pe_IP;xx];
        iter_hist_IP = [iter_hist_IP,out.iterations];
    end
    t_IP = toc
    
    [fval_IP,oi] = min(fval_hist_pe_IP);
    opt_xx_pe_IP = params_hist_pe_IP(oi,:);
    iter_IP = iter_hist_IP(oi);
    % iter_IP = sum(iter_hist_IP);
    
    
    %% fmincon IP hess optimize
    func = @(x) mix_logistic_nl_full(x,DATA,Dx,init,2,2);
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
            ,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol,'SpecifyObjectiveGradient',true);
    
    fval_hist_pe_IP_hess   = [];
    params_hist_pe_IP_hess = [];
    iter_hist_IP_hess      = [];
    tic
    for j = 1:num_opt
        [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],Aeq,beq,lb,ub,[],options1); 
        fval_hist_pe_IP_hess = [fval_hist_pe_IP_hess,ff];
        params_hist_pe_IP_hess = [params_hist_pe_IP_hess;xx];
        iter_hist_IP_hess    = [iter_hist_IP_hess,out.iterations];
    end
    t_IP_hess = toc
    
    [fval_IP_hess,oi] = min(fval_hist_pe_IP_hess);
    opt_xx_pe_IP_hess = params_hist_pe_IP_hess(oi,:);
    iter_IP_hess      = iter_hist_IP_hess(oi);
    % iter_IP_hess = sum(iter_hist_IP_hess);
    
    
    
    
    %% fmincon SQP optimize
    
    func = @(x) mix_logistic_nl_full(x,DATA,Dx,init,2,0);
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
            ,'algorithm','sqp','OptimalityTolerance',opt_tol,'StepTolerance',step_tol);
    
    fval_hist_sqp   = [];
    params_hist_sqp = [];
    iter_hist_sqp       = [];
    tic
    for j = 1:num_opt
        [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],Aeq,beq,lb,ub,[],options1); 
        fval_hist_sqp = [fval_hist_sqp,ff];
        params_hist_sqp = [params_hist_sqp;xx];
        iter_hist_sqp = [iter_hist_sqp,out.iterations];
       
    end
    t_sqp = toc
    [fval_sqp,oi] = min(fval_hist_sqp);
    opt_xx_pe_sqp = params_hist_sqp(oi,:);
    iter_sqp      = iter_hist_sqp(oi);
    % iter_sqp = sum(iter_hist_sqp);
    
    
    %% fmincon SQP grad optimize
    
    func = @(x) mix_logistic_nl_full(x,DATA,Dx,init,2,2);
    
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
            ,'algorithm','sqp','OptimalityTolerance',opt_tol,'StepTolerance',step_tol,'SpecifyObjectiveGradient',true);
    
    fval_hist_sqp_grad   = [];
    params_hist_sqp_grad = [];
    iter_hist_sqp_grad       = [];
    tic
    for j = 1:num_opt
        [xx,ff,~,out,~,g,~]  = fmincon(func,x_init(j,:),[],[],Aeq,beq,lb,ub,[],options1); 
        fval_hist_sqp_grad = [fval_hist_sqp_grad,ff];
        params_hist_sqp_grad = [params_hist_sqp_grad;xx];
        iter_hist_sqp_grad = [iter_hist_sqp_grad,out.iterations];
       
    end
    t_sqp_grad = toc
    [fval_sqp_grad,oi] = min(fval_hist_sqp_grad);
    opt_xx_pe_sqp_grad = params_hist_sqp_grad(oi,:);
    iter_sqp_grad      = iter_hist_sqp_grad(oi);
    % iter_sqp_grad = sum(iter_hist_sqp_grad);


%% Cubic

    val_cubic_f=[];
    nor_dis_i=[];   
    nor_dis_i_2=[];
    x_list_cubic_f=[];
    iter_hist_cubic = [];
    idx=[];




    tic
        for ii=1:num_opt
            % [param,val,nor_dis,nor_dis2,iter]=mysolve_mix_logistic(Dx,x_init(ii,:)',DATA,init,lb',ub',max_iter,opt_tol,step_tol);
            [param,val,nor_dis,nor_dis2,iter]=mysolve_nl(Dx,x_init(ii,:)',DATA,init,Aeq,beq,lb',ub',max_iter,opt_tol,step_tol);
            % val
            % pause
            val_cubic_f=[val_cubic_f,val(end)];
            x_list_cubic_f=[x_list_cubic_f,param];
            iter_hist_cubic = [iter_hist_cubic,iter];
            idx=[idx,ii];
        end

    t_cubic=toc

    [fval_cubic,oi_cubic] = min(val_cubic_f);
    opt_xx_pe_cubic = x_list_cubic_f(:,oi_cubic);
    iter_cubic = iter_hist_cubic(oi_cubic);
    % iter_cubic = sum(iter_hist_cubic);


    %% Record all the data

    cubic_val_hist = [cubic_val_hist;fval_cubic];
    cubic_err_hist = [cubic_err_hist;abs(opt_xx_pe_cubic-theta)./theta];
    cubic_t_hist   = [cubic_t_hist;t_cubic];
    cubic_iter_hist = [cubic_iter_hist;iter_cubic];
    fmin_IP_val_hist = [fmin_IP_val_hist;fval_IP];
    fmin_IP_err_hist = [fmin_IP_err_hist;abs(opt_xx_pe_IP-theta)./theta];
    fmin_IP_t_hist   = [fmin_IP_t_hist;t_IP];
    fmin_IP_iter_hist = [fmin_IP_iter_hist;iter_IP];
    fmin_IP_hess_val_hist = [fmin_IP_hess_val_hist;fval_IP_hess];
    fmin_IP_hess_err_hist = [fmin_IP_hess_err_hist;abs(opt_xx_pe_IP_hess-theta)./theta];
    fmin_IP_hess_t_hist   = [fmin_IP_hess_t_hist;t_IP_hess];
    fmin_IP_hess_iter_hist = [fmin_IP_hess_iter_hist;iter_IP_hess];
    fmin_sqp_val_hist = [fmin_sqp_val_hist;fval_sqp];
    fmin_sqp_err_hist = [fmin_sqp_err_hist;abs(opt_xx_pe_sqp-theta)./theta];
    fmin_sqp_t_hist   = [fmin_sqp_t_hist;t_sqp];
    fmin_sqp_iter_hist = [fmin_sqp_iter_hist;iter_sqp];
    fmin_sqp_grad_val_hist = [fmin_sqp_grad_val_hist;fval_sqp_grad];
    fmin_sqp_grad_err_hist = [fmin_sqp_grad_err_hist;abs(opt_xx_pe_sqp_grad-theta)./theta];
    fmin_sqp_grad_t_hist   = [fmin_sqp_grad_t_hist;t_sqp_grad];
    fmin_sqp_grad_iter_hist = [fmin_sqp_grad_iter_hist;iter_sqp_grad];



end





% %% Obtain the grad
% 
% 
% syms f(x,D,L,n,nx0)
% 
% f(x,D,L,n,nx0) = (L/(1+exp(-n*x+nx0)) - D)^2;
% 
% fgrad = jacobian(f(x,D,L,n,nx0),[L,n,nx0]);
% fhess = jacobian(fgrad,[L,n,nx0]);









