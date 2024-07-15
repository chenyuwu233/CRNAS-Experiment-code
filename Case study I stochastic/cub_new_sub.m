%% Cubic regurized sub problem

function [s]=cub_new_sub(g,B,sigma)
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
    if g.'*s+0.5*s.'*B*s+sigma*norm(s)^3/3>1e-3
        % keyboard
    end
end