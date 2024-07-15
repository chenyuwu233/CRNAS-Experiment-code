function hess = logistic_Hess(theta,D,x)
    hess = zeros(3,3);
    L = theta(1);
    n = theta(2);
    nx0 = theta(3);
    hess(1,1) = 2/(exp(nx0 - n*x) + 1)^2;
    hess(1,2) = -(2*x*exp(nx0 - n*x)*(D - 2*L + D*exp(nx0 - n*x)))/(exp(nx0 - n*x) + 1)^3;
    hess(1,3) = (2*exp(nx0 - n*x)*(D - L/(exp(nx0 - n*x) + 1)))/(exp(nx0 - n*x) + 1)^2 - (2*L*exp(nx0 - n*x))/(exp(nx0 - n*x) + 1)^3;
    hess(2,2) = (2*L^2*x^2*exp(2*nx0 - 2*n*x))/(exp(nx0 - n*x) + 1)^4 - (4*L*x^2*exp(2*nx0 - 2*n*x)*(D - L/(exp(nx0 - n*x) + 1)))/(exp(nx0 - n*x) + 1)^3 + (2*L*x^2*exp(nx0 - n*x)*(D - L/(exp(nx0 - n*x) + 1)))/(exp(nx0 - n*x) + 1)^2;
    hess(2,3) = (4*L*x*exp(2*nx0 - 2*n*x)*(D - L/(exp(nx0 - n*x) + 1)))/(exp(nx0 - n*x) + 1)^3 - (2*L*x*exp(nx0 - n*x)*(D - L/(exp(nx0 - n*x) + 1)))/(exp(nx0 - n*x) + 1)^2 - (2*L^2*x*exp(2*nx0 - 2*n*x))/(exp(nx0 - n*x) + 1)^4;
    hess(3,3) = (2*L^2*exp(2*nx0 - 2*n*x))/(exp(nx0 - n*x) + 1)^4 - (4*L*exp(2*nx0 - 2*n*x)*(D - L/(exp(nx0 - n*x) + 1)))/(exp(nx0 - n*x) + 1)^3 + (2*L*exp(nx0 - n*x)*(D - L/(exp(nx0 - n*x) + 1)))/(exp(nx0 - n*x) + 1)^2;
    hess(2,1) = hess(1,2);
    hess(3,1) = hess(1,3);
    hess(3,2) = hess(2,3);
end