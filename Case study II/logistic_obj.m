function ret = logistic_obj(x,DATA,Dx)
    ret = 0;
    L = x(1);
    n = x(2);
    nx0 = x(3);
    for i = 1:length(DATA)
        ret = ret + (DATA(i) - logistic(Dx(i),L,n,nx0))^2;
    end
end