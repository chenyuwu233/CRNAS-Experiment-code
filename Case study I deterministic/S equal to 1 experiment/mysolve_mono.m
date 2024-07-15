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


function [param, param_, val,nor_dis,nor_dis2,jj]=mysolve_mono(Conc,Time,NR,param,Observations_ave,MEAN_initial,num_sub,param_acc,param__acc,lb,ub,max_iter,opt_tol,step_tol)
    step_size=3e-1;
    thro=step_size^2;
    low_n=0.2;
    % lb=zeros(num_sub*5,1);     %order (p_1;p_2;alpha1;beta1;b1,...,c)
    % lb((num_sub+4):4:end)=low_n;
    % ub=Inf(num_sub*5,1);
    % ub((num_sub+2):4:end)=1;
    set_noninf=(lb~=ub);
    % set_noninf(1)=[];
    val=[];
    nor_dis=[];
    nor_dis2=[];
    adj=blkdiag(speye(4*num_sub));
    x=[vec(param)];
    t=1;
    sigma=1e4;
    sigma0=1e0;
    count=0;
    for jj=0:max_iter % Max iter
        
    [~,g_param,hess]=grad(Conc,Time,NR,1,param,Observations_ave,MEAN_initial,num_sub,2);
    hess = hess(2:end,2:end);
    if mod(jj,10)==0
        f=fun(Conc,Time,NR,1,param,Observations_ave,MEAN_initial,num_sub);
        param_=param;
        %norm(param)
        param_(3,:)=param(3,:).^(1./param(4,:));
        nor=min(norm(vec(param_-param_acc)),norm(vec(flip(param_,2)-param_acc)));
        nor2=min(norm(vec(param-param__acc)),norm(vec(flip(param,2)-param__acc)));
        if mod(jj,100)==0
        end
        val=[val,f];
        nor_dis=[nor_dis,nor];
        nor_dis2=[nor_dis2,nor2];
    end  
    ve1=min(x-lb,ub-x);
    ve=1./(ve1.^2);
    mat=blkdiag(diag(ve));
    B=adj.'*hess*adj;
    g=adj.'*[vec(g_param)];


    if norm(g)^2<opt_tol^2 % OptimalityTolerance
        1
        jj
        break
    end
    
    try
        mat_invhalf=blkdiag(diag(ve1));
    catch
        keyboard
    end
    % inv(chol(diag(ve(1:num_sub-1))+ve(num_sub)*ones(num_sub-1)))
    % mat_invhalf.'
    % g
    % keyboard
    % mat_invhalf.'*B*mat_invhalf
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
        sol_x_1=mynonconvexqp(adj.'*hess*adj,adj.'*[vec(g_param)],mat,thro,mat_invhalf);
        sol_x=sol_x_1;
        delta_x=adj*sol_x;
    end
    if sol_x.'*mat*sol_x>thro+1e-1

        % keyboard
    end
    x1=x+t*delta_x;
    ff=fun(Conc,Time,NR,1,param,Observations_ave,MEAN_initial,num_sub);
    param1=reshape(x1,4,[]);
    f2=fun(Conc,Time,NR,1,param1,Observations_ave,MEAN_initial,num_sub);

    if any(x1<0) || any (param1(2,:)>1)

        break
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
            % Pop1=x1(1:num_sub);
            param1=reshape(x1,4,[]);
            f2=fun(Conc,Time,NR,1,param1,Observations_ave,MEAN_initial,num_sub);
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
        2
        jj
        break
    end
    
    % Pop=x(1:num_sub);
    param=reshape(x,4,[]);
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