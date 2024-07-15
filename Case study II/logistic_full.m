function [obj,grad,Hess] = logistic_full(x,DATA,Dx)
    obj = 0;
    grad = zeros(3,1);
    Hess = zeros(3,3);
    for i = 1:length(DATA)
        obj = obj + (DATA(i) - logistic(Dx(i),x(1),x(2),x(3)))^2;
        grad = grad + logistic_grad(x,DATA(i),Dx(i));
        Hess = Hess + logistic_Hess(x,DATA(i),Dx(i));
    end
end