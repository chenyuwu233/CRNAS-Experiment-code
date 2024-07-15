
function Hess = mix_logistic_nl_Hess(theta, cv, DATA, init)
    
    % c = theta(end);
    % theta(end) = [];
    Theta = reshape(theta,3,[]);
    num_sub = size(Theta,2);
    Hess = zeros(3*num_sub,3*num_sub);
    prx = 0;
    grad_sub = [];
    for i = 1:num_sub
        prx = prx + init*Theta(1,i)*logistic(cv,1,Theta(2,i),Theta(3,i));
        grad_sub = [grad_sub,init*logistic_grad_sub(Theta(1,i),1,Theta(2,i),Theta(3,i),cv)];
    end
    prx = 2*(prx - DATA);
    for i = 1:num_sub
        L = 1;
        n = Theta(2,i);
        nx0 = Theta(3,i);
        p = Theta(1,i);
        exp_temp = exp(-cv*n+nx0);
        % Hess(4*i-3,4*i-3) = 2 * grad_p^2;
        % Hess(4*i-2,4*i-2) = 2 * grad_L^2;
        % Hess(4*i-1,4*i-1) = 2 * grad_n^2 + ...
        %     prx*init*p*L*cv*((-cv*exp_temp*(1+exp_temp) + 2*cv*exp_temp)/(1+exp_temp)^3);
        % Hess(4*i,4*i) = 2 * grad_nx0^2 + ...
        %     prx*init*p*(-L)*((exp_temp*(1+exp_temp) - 2*exp_temp)/(1+exp_temp)^3);
        
        
        hess_temp = prx*init*logistic_hess_sub(p,L,n,nx0,cv) + 2*grad_sub(:,i)*grad_sub(:,i)';
        hess_temp(2,:) = [];
        hess_temp(:,2) = [];
        Hess(3*i-2:3*i,3*i-2:3*i) = hess_temp;
        % logistic_hess_sub(p,L,n,nx0,cv)
        for j = i+1:num_sub
            hess_temp_ij = 2*grad_sub(:,i)*grad_sub(:,j)';
            hess_temp_ij(2,:) = [];
            hess_temp_ij(:,2) = [];
            hess_temp_ji = 2*grad_sub(:,j)*grad_sub(:,i)';
            hess_temp_ji(2,:) = [];
            hess_temp_ji(:,2) = [];
            Hess(3*i-2:3*i,3*j-2:3*j) = hess_temp_ij;
            Hess(3*j-2:3*j,3*i-2:3*i) = hess_temp_ji;
        end
    end
end