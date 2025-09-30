%% Initialization


    % number of indepdent replications of the experiment 
% 
% parpool('local',20)
% warning('off','MATLAB:integral:NonFiniteValue')

seed_num = 37
rng(seed_num)

%%
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
fmin_TR_val_hist  = [];
fmin_TR_err_hist  = [];
fmin_TR_t_hist    = [];
fmin_TR_est_hist  = [];
fmin_active_val_hist  = [];
fmin_active_err_hist  = [];
fmin_active_t_hist    = [];
fmin_active_est_hist  = [];
theta_hist = [];
x_init_hist = [];
NR = 1;
num_sub = 3;  % Change the number of subpopulation here
num_opt = 20;
num_idx = 100;
max_iter = 500;
opt_tol = 1e-6;
step_tol = 1e-6;


lb=[zeros(num_sub,1);repmat([0;0;0;0.2],num_sub,1)];
ub=[ones(num_sub,1);repmat([1;1;Inf;Inf],num_sub,1)]; % Change the upperbound here




%%

parfor i = 1:num_idx

    Time = 0:3:36;
    % DrugLevels = [0,0.0313,0.0625,0.125,0.25,0.375,0.5,1.25,2.5,3.75,5]; 
    DrugLevels = [0,logspace(-2,1,num_sub*4-1)];
    T=length(Time);
    D=length(DrugLevels);

    
    Conc = DrugLevels;
    NoiseL =0;  % lower noise level
    NoiseH =0;  % high noise level 
    
    INIT = 3290.5;
    
    
    alpha_lb = 0;alpha_ub = 0.1;
    b_lb = 0.8;b_ub = 1;
    % E_s_lb = 0.05;E_s_ub = 0.1;
    % E_m_lb = 0.1; E_m_ub = 0.5;
    % E_r_lb = 0.5;E_r_ub = 2.5;
    % E_bound = [];
    % for s = 1:num_sub
    %     E_bound = [E_bound;[mean([DrugLevels(4*s-2),DrugLevels(4*s-3)]),mean([DrugLevels(4*s),DrugLevels(4*s-1)])]];
    % end
    m_lb = 1.5;m_ub = 5;
    
    % GE_lb = [alpha_lb,b_lb,E_s_lb,m_lb,alpha_lb,b_lb,E_m_lb,m_lb,alpha_lb,b_lb,E_r_lb,m_lb];
    % GE_ub = [alpha_ub,b_ub,E_s_ub,m_ub,alpha_ub,b_ub,E_m_ub,m_ub,alpha_ub,b_ub,E_r_ub,m_ub];
    
    

    % THETA = rand(1,4*num_sub).*(GE_ub-GE_lb)+GE_lb;
    % THETA = [rand*0.3,rand*0.3,THETA,NoiseL,NoiseH];

    GE_lb = [];
    GE_ub = [];
    for s = 1:num_sub
        GE_lb = [GE_lb,[alpha_lb,b_lb,mean([DrugLevels(4*s-2),DrugLevels(4*s-3)]),m_lb]];
        GE_ub = [GE_ub,[alpha_ub,b_ub,mean([DrugLevels(4*s),DrugLevels(4*s-1)]),m_ub]];
    end

    THETA = rand(1,4*num_sub).*(GE_ub-GE_lb)+GE_lb;
    THETA = [rand(1,num_sub-1)./num_sub,THETA,NoiseL,NoiseH];
    
    
    
    
    
    Observations=generateData_param(THETA,Conc,Time,NR,num_sub,INIT);  % this function I recieved from Dr.Leder and it generates sample data 
    Observations=permute(Observations,[3,2,1]); % I permute the observations to have the data to be indexed by time, dose, replication 
    
    
    
    
    
    Observations_ave=mean(Observations,3);
    square_vari=sum((Observations-Observations_ave).^2,3);
    
    
    %% Optimization initialization
    
    low_n=0.2;
    val_fmin_IP=[];
    x_list_fmin_IP=[];
    iter_fmin_IP=[];
    val_fmin_IP_hess=[];
    x_list_fmin_IP_hess=[];
    iter_fmin_IP_hess=[];
    val_fmin_sqp=[];
    x_list_fmin_sqp=[];
    iter_fmin_sqp=[];
    val_fmin_sqp_grad=[];
    x_list_fmin_sqp_grad=[];
    iter_fmin_sqp_grad=[];
    val_fmin_TR=[];
    x_list_fmin_TR=[];
    val_fmin_active=[];
    x_list_fmin_active=[];
    
    [Pop_acc,param_acc,noise_acc] = reshape_param(THETA,num_sub);
    
    param__acc=param_acc;
    param__acc(3,:)=param_acc(3,:).^param_acc(4,:);
    
    param_ini=[];
    for ll=1:num_opt
        param=[0;0;0;0.2]+rand(4,num_sub).*[0.1;1;10;10];
        % param(3,:)=param(3,:).^param(4,:);
        Pop=ones(num_sub,1)./num_sub;
            if num_sub==2
                Pop=0.5+0.3*rand;
                Pop=[Pop;1-Pop];
            end
        param_ini=[param_ini,[Pop;vec(param)]];  % Initialization points are changed based on variable.
    end
    x_init_hist = [x_init_hist;param_ini];


        %% Fmincon IP
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
        ,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol);

    tic
    for ii=1:num_opt
        param_instance=param_ini(:,ii);
        Pop=param_instance(1:num_sub);
        param=reshape(param_instance(num_sub+1:end),4,[]);
        Aeq=[ones(1,num_sub),zeros(1,4*num_sub)];
        beq=[1];
        A=[];
        b=[];
        param_vec_0=[Pop;vec(param)];
        % fun_min=@(param_vec) mean((INIT*param_vec(1)*popfunc2(param_vec(3:6),Conc,Time.')+INIT*param_vec(2)*popfunc2(param_vec(7:10),Conc,Time.')-Observations_ave).^2,'all')/2;
        fun_min = @(theta) my_fun_obj(Conc,Time,NR,theta,Observations_ave,INIT,num_sub); % return obj gradient and hess
        [param_vec,fv,ef,output] = fmincon(fun_min,param_vec_0,[],[],Aeq,beq,lb,ub,[],options1);

        val_fmin_IP=[val_fmin_IP,fun_min(param_vec)];    % Record the objective value
        iter_fmin_IP = [iter_fmin_IP,output.iterations]; % Record the iteration number
        Pop=param_vec(1:num_sub);
        param_=reshape(param_vec(num_sub+1:end),4,[]);
        param_(3,:)=param_(3,:).^(1./param_(4,:));
        est = [Pop,param_'];
        est = sortrows(est,4)';
        % if Pop(1)>Pop(2)
        %     Pop=flip(Pop);
        %     %param=flip(param,2);
        %     param_=flip(param_,2);
        % end
        % param_vec_=[Pop;vec(param_)];

        x_list_fmin_IP=[x_list_fmin_IP,[est(1,:)';vec(est(2:end,:))]];  % Record the solution

    end
    t_IP=toc
    
    %% Fmincon IP hess
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
        ,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol,'SpecifyObjectiveGradient',true);

    tic
    for ii=1:num_opt
        param_instance=param_ini(:,ii);
        Pop=param_instance(1:num_sub);
        param=reshape(param_instance(num_sub+1:end),4,[]);
        Aeq=[ones(1,num_sub),zeros(1,4*num_sub)];
        beq=[1];
        A=[];
        b=[];
        param_vec_0=[Pop;vec(param)];
        % fun_min=@(param_vec) mean((INIT*param_vec(1)*popfunc2(param_vec(3:6),Conc,Time.')+INIT*param_vec(2)*popfunc2(param_vec(7:10),Conc,Time.')-Observations_ave).^2,'all')/2;
        fun_min = @(theta) my_fun(Conc,Time,NR,theta,Observations_ave,INIT,num_sub); % return obj gradient and hess
        [param_vec,fv,ef,output] = fmincon(fun_min,param_vec_0,[],[],Aeq,beq,lb,ub,[],options1);

        val_fmin_IP_hess=[val_fmin_IP_hess,fun_min(param_vec)];     % Record the objective
        iter_fmin_IP_hess = [iter_fmin_IP_hess, output.iterations]; % Record the iteration
        Pop=param_vec(1:num_sub);
        param_=reshape(param_vec(num_sub+1:end),4,[]);
        param_(3,:)=param_(3,:).^(1./param_(4,:));
        est = [Pop,param_'];
        est = sortrows(est,4)';
        % if Pop(1)>Pop(2)
        %     Pop=flip(Pop);
        %     %param=flip(param,2);
        %     param_=flip(param_,2);
        % end
        % param_vec_=[Pop;vec(param_)];

        x_list_fmin_IP_hess=[x_list_fmin_IP_hess,[est(1,:)';vec(est(2:end,:))]];  % Record the solution

    end
    t_IP_hess=toc

    %%  Fmincon SQP 
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
        ,'algorithm','sqp','OptimalityTolerance',opt_tol,'StepTolerance',step_tol);
    
    tic
    for ii=1:num_opt
        param_instance=param_ini(:,ii);
        Pop=param_instance(1:num_sub);
        param=reshape(param_instance(num_sub+1:end),4,[]);
        % lb=[0;0;repmat([0;0;0;0.2],2,1)];
        % ub=[1;1;repmat([1;1;Inf;10],2,1)]; % The upper bound is seted as for [ub_a,ub_b,ub_E^n,ub_n]
        Aeq=[ones(1,num_sub),zeros(1,4*num_sub)];
        beq=[1];
        A=[];
        b=[];
        param_vec_0=[Pop;vec(param)];
        % fun_min=@(param_vec) mean((INIT*param_vec(1)*popfunc2(param_vec(3:6),Conc,Time.')+INIT*param_vec(2)*popfunc2(param_vec(7:10),Conc,Time.')-Observations_ave).^2,'all')/2;
        fun_min = @(theta) my_fun_obj(Conc,Time,NR,theta,Observations_ave,INIT,num_sub); % return obj gradient and hess
        [param_vec,fv,ef,output] = fmincon(fun_min,param_vec_0,A,b,Aeq,beq,lb,ub,[],options1);
        
        val_fmin_sqp=[val_fmin_sqp,fun_min(param_vec)];    % Record the objective
        iter_fmin_sqp = [iter_fmin_sqp,output.iterations]; % Record the iteration
        Pop=param_vec(1:num_sub);
        param_=reshape(param_vec(num_sub+1:end),4,[]);
        param_(3,:)=param_(3,:).^(1./param_(4,:));

        est = [Pop,param_'];
        est = sortrows(est,4)';
        % if Pop(1)>Pop(2)
        %     Pop=flip(Pop);
        %     %param=flip(param,2);
        %     param_=flip(param_,2);
        % end
        % param_vec_=[Pop;vec(param_)];
        
        x_list_fmin_sqp=[x_list_fmin_sqp,[est(1,:)';vec(est(2:end,:))]]; % Record the solution
    
    end
    t_sqp=toc


    %%  Fmincon SQP grad
    options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',max_iter,'Display','off' ...
        ,'algorithm','sqp','SpecifyObjectiveGradient',true,'OptimalityTolerance',opt_tol,'StepTolerance',step_tol);
    
    tic
    for ii=1:num_opt
        param_instance=param_ini(:,ii);
        Pop=param_instance(1:num_sub);
        param=reshape(param_instance(num_sub+1:end),4,[]);
        % lb=[0;0;repmat([0;0;0;0.2],2,1)];
        % ub=[1;1;repmat([1;1;Inf;10],2,1)]; % The upper bound is seted as for [ub_a,ub_b,ub_E^n,ub_n]
        Aeq=[ones(1,num_sub),zeros(1,4*num_sub)];
        beq=[1];
        A=[];
        b=[];
        param_vec_0=[Pop;vec(param)];
        % fun_min=@(param_vec) mean((INIT*param_vec(1)*popfunc2(param_vec(3:6),Conc,Time.')+INIT*param_vec(2)*popfunc2(param_vec(7:10),Conc,Time.')-Observations_ave).^2,'all')/2;
        fun_min = @(theta) my_fun_grad(Conc,Time,NR,theta,Observations_ave,INIT,num_sub); % return obj gradient and hess
        [param_vec,fv,ef,output] = fmincon(fun_min,param_vec_0,A,b,Aeq,beq,lb,ub,[],options1);
        
        val_fmin_sqp_grad=[val_fmin_sqp_grad,fun_min(param_vec)];    % Record the objective
        iter_fmin_sqp_grad = [iter_fmin_sqp_grad,output.iterations]; % Record the iteration
        Pop=param_vec(1:num_sub);
        param_=reshape(param_vec(num_sub+1:end),4,[]);
        param_(3,:)=param_(3,:).^(1./param_(4,:));
                
        est = [Pop,param_'];
        est = sortrows(est,4)';
        % if Pop(1)>Pop(2)
        %     Pop=flip(Pop);
        %     %param=flip(param,2);
        %     param_=flip(param_,2);
        % end
        % param_vec_=[Pop;vec(param_)];
        
        x_list_fmin_sqp_grad=[x_list_fmin_sqp_grad,[est(1,:)';vec(est(2:end,:))]]; % Record the solution
    
    end
    t_sqp_grad=toc


    % %% Fmincon trust-region
    % options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off' ...
    %      ,'algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'StepTolerance',1e-10,'OptimalityTolerance',1e-10);
    % 
    % tic
    % for ii=1:num_opt
    %     param_instance=param_ini(:,ii);
    %     Pop=param_instance(1:num_sub);
    %     param=reshape(param_instance(num_sub+1:end),4,[]);
    %     % lb=[0;0;repmat([0;0;0;0.2],2,1)];
    %     % ub=[1;1;repmat([1;1;Inf;10],2,1)]; % The upper bound is seted as for [ub_a,ub_b,ub_E^n,ub_n]
    %     Aeq=[1,1,zeros(1,8)];
    %     beq=[1];
    %     A=[];
    %     b=[];
    %     param_vec_0=[Pop;vec(param)];
    %     % fun_min=@(param_vec) mean((INIT*param_vec(1)*popfunc2(param_vec(3:6),Conc,Time.')+INIT*param_vec(2)*popfunc2(param_vec(7:10),Conc,Time.')-Observations_ave).^2,'all')/2;
    %     fun_min = @(theta) my_fun(Conc,Time,NR,theta,Observations_ave,INIT,num_sub); % return obj gradient and hess
    %     param_vec = fmincon(fun_min,param_vec_0,A,b,Aeq,beq,lb,ub,[],options1);
    % 
    %     val_fmin_TR=[val_fmin_TR,fun_min(param_vec)];
    %     Pop=param_vec(1:2);
    %     param_=reshape(param_vec(3:end),4,[]);
    %     param_(3,:)=param_(3,:).^(1./param_(4,:));
    %     if Pop(1)>Pop(2)
    %         Pop=flip(Pop);
    %         %param=flip(param,2);
    %         param_=flip(param_,2);
    %     end
    %     param_vec_=[Pop;vec(param_)];
    % 
    %     x_list_fmin_TR=[x_list_fmin_TR,param_vec_];
    % 
    % end
    % t_TR=toc

    % %% Fmincon active set
    %     options1 = optimoptions(@fmincon,'MaxFunctionEvaluations',5990,'MaxIterations',500,'Display','off' ...
    %         ,'algorithm','active-set','StepTolerance',1e-10,'OptimalityTolerance',1e-10);
    % 
    % tic
    % for ii=1:num_opt
    %     param_instance=param_ini(:,ii);
    %     Pop=param_instance(1:num_sub);
    %     param=reshape(param_instance(num_sub+1:end),4,[]);
    %     % lb=[0;0;repmat([0;0;0;0.2],2,1)];
    %     % ub=[1;1;repmat([1;1;Inf;10],2,1)]; % The upper bound is seted as for [ub_a,ub_b,ub_E^n,ub_n]
    %     Aeq=[1,1,zeros(1,8)];
    %     beq=[1];
    %     A=[];
    %     b=[];
    %     param_vec_0=[Pop;vec(param)];
    %     fun_min=@(param_vec) mean((INIT*param_vec(1)*popfunc2(param_vec(3:6),Conc,Time.')+INIT*param_vec(2)*popfunc2(param_vec(7:10),Conc,Time.')-Observations_ave).^2,'all')/2;
    %     % fun_min = @(theta) my_fun(Conc,Time,NR,theta,Observations_ave,INIT,num_sub); % return obj gradient and hess
    %     param_vec = fmincon(fun_min,param_vec_0,A,b,Aeq,beq,lb,ub,[],options1);
    % 
    %     val_fmin_active=[val_fmin_active,fun_min(param_vec)];
    %     Pop=param_vec(1:2);
    %     param_=reshape(param_vec(3:end),4,[]);
    %     param_(3,:)=param_(3,:).^(1./param_(4,:));
    %     if Pop(1)>Pop(2)
    %         Pop=flip(Pop);
    %         %param=flip(param,2);
    %         param_=flip(param_,2);
    %     end
    %     param_vec_=[Pop;vec(param_)];
    % 
    %     x_list_fmin_active=[x_list_fmin_active,param_vec_];
    % 
    % end
    % t_active=toc


    %% Cubic reg
    
    val_i=[];
    nor_dis_i=[];   
    nor_dis_i_2=[];
    x_list=[];
    iter_cubic = [];
    best=Inf;    
    
    
    
    
    tic
    for ii=1:num_opt
        
            param_instance=param_ini(:,ii);
        Pop=param_instance(1:num_sub);
        param=reshape(param_instance(num_sub+1:end),4,[]);
        [param,param_, Pop,val,nor_dis,nor_dis2,iter]=mysolve(Conc,Time,NR,Pop,param,Observations_ave,INIT,num_sub,param_acc,param__acc,lb,ub,max_iter,opt_tol,step_tol);
%         if val(end)<best
%             best=val(end);
%             val0=val;
%         end
        val_i=[val_i,val(end)];                  % Record the objective
        iter_cubic = [iter_cubic,iter];          % Record the iteration
        nor_dis_i=[nor_dis_i,nor_dis(end)];
        nor_dis_i_2=[nor_dis_i_2,nor_dis2(end)];

        est = [Pop,param_'];
        est = sortrows(est,4)';
        
        % if Pop(1)>Pop(2)
        %     Pop=flip(Pop);
        %     % param=flip(param,2);
        %     param_=flip(param_,2);
        % end
        
        x_list=[x_list,[est(1,:)';vec(est(2:end,:))]]; % Record the solution
    end
        
            
        t_cubic=toc
    
    
    %% Recording
    
    num_best=1;
    [B,I]=mink(val_i,num_best);
    [B_fmin_IP,I_fmin_IP]=mink(val_fmin_IP,num_best);
    [B_fmin_IP_hess,I_fmin_IP_hess]=mink(val_fmin_IP_hess,num_best);
    [B_fmin_sqp,I_fmin_sqp]=mink(val_fmin_sqp,num_best);
    [B_fmin_sqp_grad,I_fmin_sqp_grad]=mink(val_fmin_sqp_grad,num_best);
    % 
    Theta_acc = [Pop_acc,reshape(param_acc,1,4*num_sub)];
    
    cubic_val_hist = [cubic_val_hist,B];
    fmin_IP_val_hist  = [fmin_IP_val_hist,B_fmin_IP];
    fmin_IP_hess_val_hist  = [fmin_IP_hess_val_hist,B_fmin_IP_hess];
    fmin_sqp_val_hist  = [fmin_sqp_val_hist,B_fmin_sqp];
    fmin_sqp_grad_val_hist  = [fmin_sqp_grad_val_hist,B_fmin_sqp_grad];
    % 
    cubic_err_hist = [cubic_err_hist,abs(x_list(:,I)-Theta_acc')./Theta_acc'];
    fmin_IP_err_hist  = [fmin_IP_err_hist,abs(x_list_fmin_IP(:,I_fmin_IP)-Theta_acc')./Theta_acc'];
    fmin_IP_hess_err_hist  = [fmin_IP_hess_err_hist,abs(x_list_fmin_IP_hess(:,I_fmin_IP_hess)-Theta_acc')./Theta_acc'];
    fmin_sqp_err_hist  = [fmin_sqp_err_hist,abs(x_list_fmin_sqp(:,I_fmin_sqp)-Theta_acc')./Theta_acc'];
    fmin_sqp_grad_err_hist  = [fmin_sqp_grad_err_hist,abs(x_list_fmin_sqp_grad(:,I_fmin_sqp_grad)-Theta_acc')./Theta_acc'];
    % 
    cubic_t_hist   = [cubic_t_hist,t_cubic];
    fmin_IP_t_hist    = [fmin_IP_t_hist,t_IP];
    fmin_IP_hess_t_hist    = [fmin_IP_hess_t_hist,t_IP_hess];
    fmin_sqp_t_hist    = [fmin_sqp_t_hist,t_sqp];
    fmin_sqp_grad_t_hist    = [fmin_sqp_grad_t_hist,t_sqp_grad];
    % 
    cubic_est_hist = [cubic_est_hist,x_list(:,I)];
    fmin_IP_est_hist = [fmin_IP_est_hist,x_list_fmin_IP(:,I_fmin_IP)];
    fmin_IP_hess_est_hist = [fmin_IP_hess_est_hist,x_list_fmin_IP_hess(:,I_fmin_IP_hess)];
    fmin_sqp_est_hist = [fmin_sqp_est_hist,x_list_fmin_sqp(:,I_fmin_sqp)];
    fmin_sqp_grad_est_hist = [fmin_sqp_grad_est_hist,x_list_fmin_sqp_grad(:,I_fmin_sqp_grad)];
    %
    cubic_iter_hist = [cubic_iter_hist,iter_cubic(I)];
    fmin_IP_iter_hist = [fmin_IP_iter_hist,iter_fmin_IP(I_fmin_IP)];
    fmin_IP_hess_iter_hist = [fmin_IP_hess_iter_hist,iter_fmin_IP_hess(I_fmin_IP_hess)];
    fmin_sqp_iter_hist = [fmin_sqp_iter_hist,iter_fmin_sqp(I_fmin_sqp)];
    fmin_sqp_grad_iter_hist = [fmin_sqp_grad_iter_hist,iter_fmin_sqp_grad(I_fmin_sqp_grad)];
    %
    cubic_number_guess_hist = [cubic_number_guess_hist,sum(val_i<1)];
    fmin_IP_number_guess_hist = [fmin_IP_number_guess_hist,sum(val_fmin_IP<1)];
    fmin_IP_hess_number_guess_hist = [fmin_IP_hess_number_guess_hist,sum(val_fmin_IP_hess<1)];
    fmin_sqp_number_guess_hist = [fmin_sqp_number_guess_hist,sum(val_fmin_sqp<1)];
    fmin_sqp_grad_number_guess_hist = [fmin_sqp_grad_number_guess_hist,sum(val_fmin_sqp_grad<1)];

    theta_hist = [theta_hist,Theta_acc'];

end

name = strcat('Result/opt',num2str(num_opt),'_',num2str(num_idx),'_uc_ng_sub',num2str(num_sub),'_no_rounding.mat');
save(name)