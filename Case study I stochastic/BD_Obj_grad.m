% This is an objective function with Branching process, for 2 subpopulation
%
% 
% Input: 
%   - Params: (1 x 13 vector) all necessary parameters ([p_i,r_i, d_i, b_i, E_i^n_i, n_i], c)
%   - DATA: (NT x ND x NR matrix) observation
%   - Conc: (1 x ND vector) concentration levels
%   - Time: (1 x NT vector) of time points
%   - n: (scalar) initial total cell number.
%   - cmd: (string) command to obtain the desired results
%
% Output:
%   - obj: (scalar) objective value
%   - grad: (1 x 6*S+1 vector) gradient at Params
%   - Hess: (6*S+1 x 6*S+1 matrix) Hessian matrix at Params
%
% Requirment:
%   - mu(d,t,theta)
%   - sig(d,t,theta)


function [obj,grad] = BD_Obj_grad(Params,DATA,Conc,Time,n,cmd)    
    %% Separate the parameter for each sub-group 
    pl   = length(Params);
    c      = Params(end);
    Params(end) = [];
    PMAT   = reshape(Params,6,[])'; % S x 6 parameters matrix
    num_sub = size(PMAT,1);

    
    
    %% Record the dimension
    
    [NT,ND,NR] = size(DATA);
    
%     %% Record the mean and variance matrix: NT x ND
%     
%     [mean_mat,var_mat] = BD_MS_vec(PMAT,Time,Conc,n.*p_vec);
%     mean_mat = mean_mat(2:end,:);
%     var_mat = var_mat(2:end,:);
%     
%     %% Record the objective
%     obj = -NR * sum(log(1./sqrt(2*pi.*(var_mat+c^2))),"all");
%     for i = 1:NR
%         obj = obj + sum((DATA(2:end,:,i)-mean_mat).^2./(2*(var_mat+c^2)),"all");
%     end

    %% Record the objective, gradient, Hessian

    obj  = 0;
    grad = zeros(1,pl);
    for i = 1:length(Conc)
        for j = 2:length(Time)

            [obj_mu_s, grad_mu_s] = mu(Conc(i),Time(j),PMAT(1,:),cmd);
            grad_mu_s = squeeze(grad_mu_s)';
            [obj_sig_s, grad_sig_s] = sig(Conc(i),Time(j),PMAT(1,:),cmd);
            [obj_mu_r, grad_mu_r] = mu(Conc(i),Time(j),PMAT(2,:),cmd);
            grad_mu_r = squeeze(grad_mu_r)';
            [obj_sig_r, grad_sig_r] = sig(Conc(i),Time(j),PMAT(2,:),cmd);

            TV = n*(obj_sig_s + obj_sig_r) + c^2;
            TM = n*(obj_mu_s + obj_mu_r);

            %% Record Objective
            obj = obj - NR*log(1/sqrt(2*pi*(TV)));
            obj = obj + sum((DATA(j,i,:)-TM).^2./(2*TV),'all');


            %% Formulate pop Gradient

            grad_mu = [grad_mu_s,grad_mu_r];
            grad_sig = [grad_sig_s,grad_sig_r];

            %% Record Gradient

            grad(1:12) = grad(1:12) + n*NR/(2*TV).*grad_sig;
            grad(13) = grad(13) + NR*c/TV;

            for k = 1:NR
                grad(1:12) = grad(1:12) - n*(DATA(j,i,k) - TM)/TV.*grad_mu;
                grad(1:12) = grad(1:12) - n*(DATA(j,i,k) - TM)^2/(2*TV^2).*grad_sig;
                grad(13)  = grad(13) - c*(DATA(j,i,k) - TM)^2/(TV^2);
            end
        end
    end

end