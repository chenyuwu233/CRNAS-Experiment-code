%% Non convex quadratic programing

function [x]=mynonconvexqp(A,b,Q,c,R_inv)
    x1=zeros(size(A,2),1);
    set_noninf=all(~isinf(Q),2);
    A=A(set_noninf,set_noninf);
    b=b(set_noninf);
    Q=Q(set_noninf,set_noninf);
    R_inv=R_inv(set_noninf,set_noninf);
    
    if nargin <5
        [R,flag]=chol(Q);
        if flag~=0
            [R,flag]=chol(Q+1e-4*eye(size(Q)));
        end
        R_inv=inv(R);
    end
    [V,d_] = eig(R_inv.'*A*R_inv, 'vector');
    [d_, ind] = sort(d_,'ascend');
    V = V(:, ind);
    if d_(1)>0
        x=-A\b;
        if x.'*Q*x <= c
            x1(set_noninf)=x;
          
            x=x1;
            return 
        end
    end
    b_=V.'*R_inv.'*b;
    num=sum(d_==d_(1));
    c1=sum((b_((num+1):end)./(d_((num+1):end)-d_(1))).^2);
    c2=b_(1:num).'*b_(1:num);
    c3=b_.'*b_;
    if c2>1e-8 || c1>c+1e-6
        mu_l=max(sqrt(c2/c)-d_(1),0);
        mu_r=sqrt(c3/c)-d_(1);
        
        for jj=1:40
            mu=(mu_l+mu_r)/2;
            x_=-b_./(d_+mu);
            f=x_.'*x_-c;
            if f>0
                mu_l=mu;
            else
                mu_r=mu;
            end
            
            if abs(f)<1e-7
                break
            end
        end
        %mu=(mu_l+mu_r)/2;
        x=-(A+mu*Q)\b;
    else
        x_=zeros(size(b_));
        x_((num+1):end)=-b_((num+1):end)./(d_((num+1):end)-d_(1));
        x_(1)=sqrt(c-x_.'*x_);
        x=R_inv*V*x_;
        mu=-d_(1);
    end
    % if mu<0
    %     keyboard
    % end
    % if abs(x.'*Q*x-c)>1e-6
    %     keyboard
    % end
    % if norm((A+mu*Q)*x+b)>1e-5
    %     keyboard
    % end
    % if any(eig(A+mu*Q)<-1e-5)
    %     keyboard
    % end
    x1(set_noninf)=x;
    x=x1;
end