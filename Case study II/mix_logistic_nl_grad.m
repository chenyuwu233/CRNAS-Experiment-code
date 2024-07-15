function grad = mix_logistic_nl_grad(theta, cv, DATA, init)
    
    % c = theta(end);
    % theta(end) = [];
    Theta = reshape(theta,3,[]);
    num_sub = size(Theta,2);
    grad = zeros(3*num_sub,1);
    prx = 0;
    for i = 1:num_sub
        prx = prx + init*Theta(1,i)*logistic(cv,1,Theta(2,i),Theta(3,i));
    end
    prx = 2*(prx - DATA);
    for i = 1:num_sub
        L = 1;
        n = Theta(2,i);
        nx0 = Theta(3,i);
        p = Theta(1,i);
        % grad(4*i-3) = prx*init*logistic(cv,L,n,nx0);
        % grad(4*i-2) = prx*init*p/(1+exp(-cv*n+nx0));
        % grad(4*i-1) = prx*init*p*(L*cv*exp(-cv*n+nx0)/(1+exp(-cv*n+nx0))^2);
        % grad(4*i)   = prx*init*p*(-L*exp(-cv*n+nx0)/(1+exp(-cv*n+nx0))^2);
        grad_temp = logistic_grad_sub(p,L,n,nx0,cv);
        grad_temp(2) = [];
        grad(3*i-2:3*i) = prx*init*grad_temp;
    end
end