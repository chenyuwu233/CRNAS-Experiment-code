%% expected population from parameter vector b

function [pop] = popfunc2(b,x, T)

pop = exp(T* (b(1)+log(b(2)+ (1-b(2))./(1+x.^b(4)./b(3)))));  % x.^b(4).*b(3)
 
end