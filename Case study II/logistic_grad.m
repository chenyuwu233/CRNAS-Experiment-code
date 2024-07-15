function grad = logistic_grad(theta,D,x)
    grad = zeros(3,1);
    L = theta(1);
    n = theta(2);
    nx0 = theta(3);
    grad(1) = -(2*(D - L + D*exp(nx0 - n*x)))/(exp(nx0 - n*x) + 1)^2;
    grad(2) = -(2*L*x*exp(nx0 - n*x)*(D - L + D*exp(nx0 - n*x)))/(exp(nx0 - n*x) + 1)^3;
    grad(3) = (2*L*exp(nx0 - n*x)*(D - L + D*exp(nx0 - n*x)))/(exp(nx0 - n*x) + 1)^3;
end