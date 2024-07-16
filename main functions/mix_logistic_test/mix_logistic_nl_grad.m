function grad = mix_logistic_nl_grad(theta,  DATA, cv, inputs)
    init = inputs.init;
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
        grad_temp = logistic_grad_sub(p,L,n,nx0,cv);
        grad_temp(2) = [];
        grad(3*i-2:3*i) = prx*init*grad_temp;
    end
end