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


function [param, param_, Pop, val,nor_dis,nor_dis2]=mysolve_lbub(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub,param_acc,param__acc,lb,ub)
    step_size=3e-1;
    thro=step_size^2;
    low_n=0.2;
    % lb=zeros(num_sub*5,1); %order (p_1;p_2;alpha1;beta1;b1,...)
    % lb((num_sub+4):4:end)=low_n;
    % ub=Inf(num_sub*5,1);
    % ub((num_sub+2):4:end)=1;
    val=[];
    nor_dis=[];
    nor_dis2=[];
    adj=blkdiag([eye(num_sub-1);-ones(1,num_sub-1)],speye(4*num_sub));
    x=[Pop;vec(param)];
    n=5*num_sub;
    t=1;
    sigma=1e4;
    sigma0=1e0;
    count=0;
    
    for jj=0:100
    
        [g_p,g_param,hess]=grad(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub,2);
        if mod(jj,10)==0
            f=fun(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub);
            param_=param;
            param_(3,:)=param(3,:).^(1./param(4,:));
            nor=min(norm(vec(param_-param_acc)),norm(vec(flip(param_,2)-param_acc)));
            nor2=min(norm(vec(param-param__acc)),norm(vec(flip(param,2)-param__acc)));
            if mod(jj,100)==0
                %fprintf('%d %.2e %.2e %.2e %.2e %.2e\n', jj, 1-sum(Pop),norm(g_param),f,nor,nor2);
            end
            %1-sum(Pop)
            val=[val,f];
            nor_dis=[nor_dis,nor];
            nor_dis2=[nor_dis2,nor2];
        end
        % 
        % ve=[1./param(1,:).^2;1./min(param(2,:),1-param(2,:)).^2;1./param(3,:).^2;1./(param(4,:)-low_n).^2];
        % mat=blkdiag(diag(1./Pop(1:end-1).^2)+1/Pop(end)^2*ones(num_sub-1),diag(vec(ve)));
        % adj=blkdiag([eye(num_sub-1);-ones(1,num_sub-1)],speye(4*num_sub));
        % 
        ve1=min(x-lb,ub-x);
        ve=1./(ve1.^2);
        mat=blkdiag(diag(ve(1:num_sub-1))+ve(num_sub)*ones(num_sub-1),diag(ve((num_sub+1):end)));
        mat1=mat;
        mat1(isinf(mat1))=0;
        B=adj.'*hess*adj;
        g=adj.'*[g_p;vec(g_param)];
        %mat_invhalf=blkdiag(inv(chol(diag(1./Pop(1:end-1).^2)+1/Pop(end)^2*ones(num_sub-1))),diag(vec(1./sqrt(ve))));
        mat_invhalf=blkdiag(inv(chol(diag(ve(1:num_sub-1))+ve(num_sub)*ones(num_sub-1))),diag(ve1((num_sub+1):end)));
        
        s=mat_invhalf*cub_new_sub(mat_invhalf.'*g,mat_invhalf.'*B*mat_invhalf,sigma);
        delta_x=adj*s;
        if s.'*mat1*s<thro
            delta_x=adj*s;
            count=count+1;
            sol_x=s;
        else
            sol_x_1=mynonconvexqp(adj.'*hess*adj,adj.'*[g_p;vec(g_param)],mat,thro,mat_invhalf);
            sol_x=sol_x_1;
            delta_x=adj*sol_x;
        end
        x1=x+t*delta_x;
        ff=fun(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub);
        Pop1=x1(1:num_sub);
        param1=reshape(x1(num_sub+1:end),4,[]);
        f2=fun(Conc,Time,NR,Pop1,param1,Observations_ave,MEAN_initial,num_sub);
        if any(x1<0) || any (param1(2,:)>1)
            keyboard
            error('Wrong 1')
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
                if s.'*mat1*s<thro
                    delta_x=adj*s;
                    count=count+1;
                    sol_x=s;
                else
                
                    sol_x=sol_x_1;
                    delta_x=adj*sol_x;
                end
                x1=x+t*delta_x;
                Pop1=x1(1:num_sub);
                param1=reshape(x1(num_sub+1:end),4,[]);
                f2=fun(Conc,Time,NR,Pop1,param1,Observations_ave,MEAN_initial,num_sub);
            else
                break
            end
            if jjj==20 && f2>ff+1e-7
                keyboard
            end
        end
        sigma=max(sigma/2,sigma0);
        %x=mylinesearch(x,delta_x,Conc,Time,NR,Observations_ave,MEAN_initial,num_sub);
        x=x+t*delta_x;
        
        Pop=x(1:num_sub);
        param=reshape(x(num_sub+1:end),4,[]);
    end
    %sigma
    %count
end