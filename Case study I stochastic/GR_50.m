function ret = GR_50(a,b,E,n,d)
    r_m = (a + a + log(b + (1 - b)/(1 + (d/E)^n)))/2;
    ret = E * ((exp(r_m - a) - 1)/(b - exp(r_m - a)))^(1/n);
end