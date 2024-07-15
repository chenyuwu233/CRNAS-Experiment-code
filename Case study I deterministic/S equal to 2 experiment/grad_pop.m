%%  Gets the 0-2 order information from model using Hill growth rate function parametrization

function [pop,g1,g2,g3,g4,hess] = grad_pop(b,x, T,opt)
%%  Gets population size from model using Hill growth rate function parametrization
%b = parameters vector
%b(1) alpha
%b(2) b
%b(3) E
%b(4) n

%x = drug concentration
%T = time

log_x=mylog(x);
%x concentration
if nargin <4
de=x.^b(4)/b(3);
de2=b(2)+ (1-b(2))./(1+de);
pop = exp(T.* (b(1)+log(de2)));
g1=T.*pop;
g2=-T.*pop.*(1./(de + 1) - 1)./de2;
g3=-(T.*b(3).*de.*pop.*(b(2) - 1)) ./(b(3)^2.*(de + 1).^2.*de2);
g4=(T.*de.*pop.*log_x.*(b(2) - 1))./((de + 1).^2.*de2);
else
    de=x.^b(4)./b(3);
    de2=b(2)+ (1-b(2))./(1+de);
    pop = exp(T.* (b(1)+log(de2)));
    g1=T.*pop;
    temp1=(1./(de + 1) - 1)./de2;
    g2=-T.*pop.*temp1;
    temp2=de./((de + 1).^2.*de2);
    temp22=temp2.^2;
    temp3=T.*temp2.*pop;
    g3=-(b(2) - 1)/b(3)*temp3;
    temp4=temp3.*log_x;
    g4=(b(2) - 1)*temp4;

    temp5=(T.^2-T).*temp22.*pop;
    temp6=temp3./(de+1);
    hess11=T.*g1;
    hess22=(T.^2-T).*pop.*temp1.^2;
    hess33=2*(b(2) - 1)/b(3)^2*temp6 + (b(2) - 1)^2/(b(3)^2)*temp5;
    hess44=(((b(2) - 1)*temp6.*(1-de)) + (b(2) - 1).^2*temp5).*log_x.^2;
    hess12=T.*g2;
    hess13=T.*g3;
    hess14=T.*g4;
    hess23=((T-1).*temp3.*temp1.*(b(2) - 1))./(b(3)) - 1/b(3)*temp3;
    hess24=temp4 - ((T-1).*temp4.*temp1.*(b(2) - 1));
    hess34=(b(2) - 1)/b(3)*temp6.*(de-1).*log_x - (b(2) - 1)^2/b(3)*(temp5.*log_x);

    hess=cat(4,cat(3,hess11,hess12,hess13,hess14),...
    cat(3,hess12,hess22,hess23,hess24),...
    cat(3,hess13,hess23,hess33,hess34),...
    cat(3,hess14,hess24,hess34,hess44));



end
end