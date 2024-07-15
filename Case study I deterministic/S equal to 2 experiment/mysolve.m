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


function [param, param_, Pop, val,nor_dis,nor_dis2,jj]=mysolve(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub,param_acc,param__acc,lb,ub,max_iter,opt_tol,step_tol)
    step_size=3e-1;
    thro=step_size^2;
    low_n=0.2;
    % lb=zeros(num_sub*5,1);     %order (p_1;p_2;alpha1;beta1;b1,...,c)
    % lb((num_sub+4):4:end)=low_n;
    % ub=Inf(num_sub*5,1);
    % ub((num_sub+2):4:end)=1;
    set_noninf=(lb~=ub);
    set_noninf(1)=[];
    val=[];
    nor_dis=[];
    nor_dis2=[];
    adj=blkdiag([eye(num_sub-1);-ones(1,num_sub-1)],speye(4*num_sub));
    x=[Pop;vec(param)];
    n=5*num_sub;
    t=1;
    sigma=1e4;
    sigma0=1e-2;
    count=0;
    x_0=x;
    for jj=0:max_iter % Max iter
        % sigma
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
%         ve=[1./param(1,:).^2;1./min(param(2,:),1-param(2,:)).^2;1./param(3,:).^2;1./(param(4,:)-low_n).^2];
%         mat=blkdiag(diag(1./Pop(1:end-1).^2)+1/Pop(end)^2*ones(num_sub-1),diag(vec(ve)));
%         adj=blkdiag([eye(num_sub-1);-ones(1,num_sub-1)],speye(4*num_sub));
%         
    ve1=min(x-lb,ub-x);
    % ve1
    ve=1./(ve1.^2);
    mat=blkdiag(diag(ve(1:num_sub-1))+ve(num_sub)*ones(num_sub-1),diag(ve((num_sub+1):end)));
    

    
    % mat
    % keyboard
    
    
    B=adj.'*hess*adj;
    g=adj.'*[g_p;vec(g_param)];
    % size(B)
    % pause

%     set_noninf=all(~isinf(B),2);
%     B=B(set_noninf,set_noninf);
% g=g(set_noninf);
% Q=Q(set_noninf,set_noninf);
% R_inv=R_inv(set_noninf,set_noninf);


    %mat_invhalf=blkdiag(inv(chol(diag(1./Pop(1:end-1).^2)+1/Pop(end)^2*ones(num_sub-1))),diag(vec(1./sqrt(ve))));

    
    try
        mat_invhalf=blkdiag(inv(chol(diag(ve(1:num_sub-1))+ve(num_sub)*ones(num_sub-1))),diag(ve1((num_sub+1):end)));
    catch
        % keyboard
    end

    if norm(mat_invhalf.'*g)^2<opt_tol^2 % OptimalityTolerance
        % 1
        % jj
        break
    end
    
    s=mat_invhalf*cub_new_sub(mat_invhalf.'*g,mat_invhalf.'*B*mat_invhalf,sigma,set_noninf);
    
    
    if ~isreal(s)
        % keyboard
    end
    delta_x=adj*s;
    if s.'*mat*s<thro
        delta_x=adj*s;
        count=count+1;
        sol_x=s;
    else
        sol_x_1=mynonconvexqp(adj.'*hess*adj,adj.'*[g_p;vec(g_param)],mat,thro,mat_invhalf);
        % if length(sol_x_1)<size(mat,1)
        %     inf_idx = linspace(1,size(mat,1),size(mat,1))*(1-all(~isinf(mat),2));
        %     if inf_idx == 1
        %         sol_x_1 = [0;sol_x_1];
        %     elseif inf_idx == size(mat,1)
        %         sol_x_1 = [sol_x_1;0];
        %     else
        %         temp_sol = zeros(size(mat,1),1);
        %         temp_sol(1:inf_idx-1) = sol_x_1(1:inf_idx-1);
        %         temp_sol(inf_idx+1:end) = sol_x_1(inf_idx:end);
        %         sol_x_1 = temp_sol;
        %     end
        % end
        sol_x=sol_x_1;
        delta_x=adj*sol_x;
    end
    if sol_x.'*mat*sol_x>thro+1e-1

        % keyboard
    end
    x1=x+t*delta_x;
    ff=fun(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub);
    Pop1=x1(1:num_sub);
    param1=reshape(x1(num_sub+1:end),4,[]);
    f2=fun(Conc,Time,NR,Pop1,param1,Observations_ave,MEAN_initial,num_sub);

    if any(x1<0) || any (param1(2,:)>1)
        % keyboard
        % x1
        % param1(2,:)
        break
        % error('Wrong 1')
    end
    %%%%%%%%%%%
    %sigma1=sigma;
    for jjj=1:100
        if f2>ff+1e-7%%%%+sol_x.'*g+0.5*sol_x.'*B*sol_x + sigma*(sol_x.'*mat*sol_x)^1.5/3 %%%%%%%
            sigma=2*sigma;
%                 jjj
%                 sigma
            s=mat_invhalf*cub_new_sub(mat_invhalf.'*g,mat_invhalf.'*B*mat_invhalf,sigma,set_noninf);
            delta_x=adj*s;
            if s.'*mat*s<thro
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
        if jjj==100 && f2>ff+1e-7
            % keyboard
        end
    end
    if ~isreal(delta_x)
        % keyboard
    end
    sigma=max(sigma/2,sigma0);
    %x=mylinesearch(x,delta_x,Conc,Time,NR,Observations_ave,MEAN_initial,num_sub);
    x=x+t*delta_x;
    if norm(delta_x)<step_tol
        % 2
        % jj
        break
    end
    
    Pop=x(1:num_sub);
    param=reshape(x(num_sub+1:end),4,[]);
    end
    %sigma
        %count
end
function [s]=cub_new_sub(g,B,sigma,set_noninf)
    s1=zeros(size(g));
    B=B(set_noninf,set_noninf);
    g=g(set_noninf);


        [V,d_] = eig(B, 'vector');
        [d_, ind] = sort(d_,'ascend');
        V = V(:, ind);
        g1=V.'*g;
        num=sum(d_==d_(1));

        
        c2=g1(1:num).'*g1(1:num);
if c2<1e-10 && d_(1)<=0
    s_=zeros(size(d_));
    s_((num+1):end)=-g1((num+1):end)./(d_((num+1):end)-d_(1));
    lambda=-d_(1);
    if norm(s_)<=lambda/sigma
    s_(1)=sqrt((lambda/sigma)^2-s_.'*s_);
    s=V*s_;
        if g.'*s+0.5*s.'*B*s+sigma*norm(s)^3/3>1e-6
            % keyboard
        end
    s1(set_noninf)=s;
    s=s1;
    return
    end
end
        lambda_delta=1;
        lambda=max(0,-d_(1))+lambda_delta;
        
        for ii=1:30
            if sqrt(sum((g1./(lambda+d_)).^2))>lambda/sigma
                break
            end
            lambda_delta=lambda_delta/4;
            lambda=max(0,-d_(1))+lambda_delta;
        end
        if ii==30
            % keyboard
        end
        
        for ll=0:100
            
        %L=chol(B+lambda*eye(size(B,1))).';
        %s=-L.'\(L\g);
        %w=L\s;
        s=-V*(diag(1./(d_+lambda))*(V.'*g));
        w=diag(1./sqrt(d_+lambda))*(V.'*s);
        delta_lambda=lambda*(norm(s)-lambda/sigma)/(norm(s)+lambda/sigma*(lambda*norm(w)^2/norm(s)^2));
        
        lambda=lambda+delta_lambda;
        %delta_lambda
            if abs(norm(s)-lambda/sigma)<1e-8
                break
            end
        end
        if g.'*s+0.5*s.'*B*s+sigma*norm(s)^3/3>1
            % keyboard
        end
            s1(set_noninf)=s;
    s=s1;
end