% This is an solver function utilizing the cubic regurized algorithm.
%
% 
% Input: 
%   - Params: (1 x 2*S+1 vector) all necessary parameters ([p_i,r_i, d_i, b_i, E_i^n_i, n_i], c)
%   - DATA: (NT x ND x NR matrix) observation
%   - Conc: (1 x ND vector) concentration levels
%   - Time: (1 x NT vector) of time points
%   - MEAN_initial: (scalar) initial total cell number
%   - num_sub: (scalar) number of sub-population: S
%   - lb: (1 x 2*S+1 vector) lower bound ([p_i],[r_i,d_i,b_i,E_i^n_i,n_i],c)
%   - ub: (1 x 2*S+1 vector) upper bound ([p_i],[r_i,d_i,b_i,E_i^n_i,n_i],c)
%   - max_iter: (scalar) Maximum iteration
%
% Output:
%   - Params: (1) objective value
%   - grad: (1 x 6*S+1 vector) gradient at Params
%   - Hess: (6*S+1 x 6*S+1 matrix) Hessian matrix at Params
%
% Requirment:
%   - mu(d,t,theta)
%   - sig(d,t,theta)



function [Params, val,nor_dis,nor_dis2,jj]=mysolve_EP(Conc,Time,NR,Params,DATA,MEAN_initial,num_sub,lb,ub,max_iter,optimality_tol,step_tol)
    step_size=3e-1;
    thro=step_size^2;
    low_n=0.2;
    % lb=zeros(size(Params));     %order (p_1;p_2;alpha1;beta1;b1,...,c)
    % lb((num_sub+5):5:end)=low_n;
    % ub=Inf(size(Params));
    % ub((num_sub+1):5:end)=1;
    % ub((num_sub+2):5:end)=1;
    % ub((num_sub+3):5:end)=1;
    val=[];
    nor_dis=[];
    nor_dis2=[];
      %heristic subproblem solver
    order=reshape(1:6*num_sub,6,[]);
    order=[order(1,:),vec(order(2:end,:)).',6*num_sub+1] ;     
    perm=sparse(1:(6*num_sub+1),order,1,6*num_sub+1,6*num_sub+1) ;  
    adj=blkdiag([eye(num_sub-1);-ones(1,num_sub-1)],speye(5*num_sub+1));
    t=1;
    sigma=1e4;
    sigma0=1e0;
    count=0;
            
        for jj=0:max_iter
        [obj,grad,Hess] = BD_Obj_full(Params,DATA,Conc,Time,MEAN_initial,'Hess');
        if jj == 0
            grad_x0 = grad;
            relative_opt  = norm(grad_x0,inf);
        end
        % [obj,grad,Hess] = BD_Obj_full_finite(Params,DATA,Conc,Time,MEAN_initial,1e-10,1e-6);
        grad=perm*grad.';
        Hess=perm*Hess*perm.';
        if mod(jj,10)==0
           % obj = BD_Obj_full(Params,DATA,Conc,Time,MEAN_initial,1);   
    % param_=param;
    % param_(3,:)=param(3,:).^(1./param(4,:));
    % nor=min(norm(vec(param_-param_acc)),norm(vec(flip(param_,2)-param_acc)));
    % nor2=min(norm(vec(param-param__acc)),norm(vec(flip(param,2)-param__acc)));
    %             
            val=[val,obj];
    % nor_dis=[nor_dis,nor];
    % nor_dis2=[nor_dis2,nor2];
        end
        param=reshape(Params(1:end-1),6,[]);
        Pop=param(1,:).';
    %         ve1=[min(param(2:4,:),1-param(2:4,:));param(5,:);param(6,:)-low_n];
    %         ve=1./(ve1.^2);
    %         %ve=[1./param(2,:).^2;1./min(param(3,:),1-param(3,:)).^2;1./param(3,:).^2;1./(param(4,:)-low_n).^2];
    %         mat=blkdiag(diag(1./Pop(1:end-1).^2)+1/Pop(end)^2*ones(num_sub-1),diag(vec(ve)),1./Params(end).^2);
    %         
        
        Param_order=perm*Params;
        ve1=min(Param_order-lb,ub-Param_order);
        ve=1./(ve1.^2);
        mat=blkdiag(diag(ve(1:num_sub-1))+ve(num_sub)*ones(num_sub-1),diag(ve((num_sub+1):end)));
        B=adj.'*Hess*adj;
        if norm(B.'-B)>0
            % keyboard;
        end
        g=adj.'*grad;
        if norm(grad(1:num_sub)-mean(grad(1:num_sub)))^2+norm(grad(num_sub+1:end))^2<optimality_tol^2 || norm(grad(1:num_sub)-mean(grad(1:num_sub)))^2+norm(grad(num_sub+1:end))^2 <=(relative_opt*optimality_tol)^2% Optimality condition
            1
            jj
            break
        end


        try
        %mat_invhalf=blkdiag(inv(chol(diag(1./Pop(1:end-1).^2)+1/Pop(end)^2*ones(num_sub-1))),diag(vec(1./sqrt(ve))),Params(end));
        mat_invhalf=blkdiag(inv(chol(diag(ve(1:num_sub-1))+ve(num_sub)*ones(num_sub-1))),diag(ve1((num_sub+1):end)));
        catch
            % keyboard;
        end
        s=mat_invhalf*cub_new_sub(mat_invhalf.'*g,mat_invhalf.'*B*mat_invhalf,sigma);
        delta_x=adj*s;
        if s.'*mat*s<thro
            delta_x=adj*s;
            count=count+1;
            sol_x=s;
        else
            sol_x_1=mynonconvexqp(adj.'*Hess*adj,adj.'*grad,mat,thro,mat_invhalf);
            sol_x=sol_x_1;
            delta_x=adj*sol_x;
        end
        x1=Params+t*perm.'*delta_x;
        %ff=fun(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub);
        ff=obj;
        x1_reshape=reshape(x1(1:end-1),6,[]);
        f2= BD_Obj(x1,DATA,Conc,Time,MEAN_initial,'Obj'); %only function value is need
        
        %keyboard
        if any(perm*x1>ub) || any (perm*x1<lb)
            % keyboard
            % error('Wrong 1')
        end
        %%%%%%%%%%%
        %sigma1=sigma;
        for jjj=1:20
            if f2>ff+1e-7%%%%+sol_x.'*g+0.5*sol_x.'*B*sol_x + sigma*(sol_x.'*mat*sol_x)^1.5/3 %%%%%%%
                sigma=2*sigma;
    % jjj
    % sigma
                s=mat_invhalf*cub_new_sub(mat_invhalf.'*g,mat_invhalf.'*B*mat_invhalf,sigma);
                delta_x=adj*s;
                if s.'*mat*s<thro
                    delta_x=adj*s;
                    count=count+1;
                    sol_x=s;
                else
                   
                    sol_x=sol_x_1;
                    delta_x=adj*sol_x;
                end
                x1=Params+t*perm.'*delta_x;
                
                f2= BD_Obj(x1,DATA,Conc,Time,MEAN_initial,'Obj'); %only function value is needed
            else
                break
            end
            
        end
            if ~isreal(delta_x)
                % keyboard;
            end
        sigma=max(sigma/2,sigma0);
        %x=mylinesearch(x,delta_x,Conc,Time,NR,Observations_ave,MEAN_initial,num_sub);
        Params_pre = Params;
        Params=Params+t*perm.'*delta_x;
        if norm(perm.'*delta_x)<step_tol || norm(perm.'*delta_x)<step_tol*(1+norm(Params_pre))  % Step tol condition
            2
            jj
            break
        end
            if ~isreal(Params)
                % keyboard;
            end
        end
%sigma
    %count
end