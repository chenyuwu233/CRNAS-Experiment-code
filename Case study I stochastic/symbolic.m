vps = 0.4;
vpr = 1-vps;
vbetas = 0.56;
vnus = 0.54;
vbs = 0.87;
vEs = 0.08;
vms = 2.5;
vbetar = 0.13;
vnur = 0.08;
vbr = 0.97;
vEr = 1;
vmr = 1.5;
c = 0;
Time = 0:3:3;
% Conc = [0,0.25,0.375,0.5,1.25,2.5];  
Conc = 0.375;
init = 1000;
NT = length(Time);
ND = length(Conc);
NR = 1;
theta = [vps,vbetas,vnus,vbs,vEs,vms,vpr,vbetar,vnur,vbr,vEr,vmr,c];
PMAT = [vps,vbetas,vnus,vbs,vEs^vms,vms;
        vpr,vbetar,vnur,vbr,vEr^vmr,vmr];
DATA = sto_gen_bd(NR,Conc,Time,init,theta);
theta_c = [vps,vbetas,vnus,vbs,vEs^vms,vms,vpr,vbetar,vnur,vbr,vEr^vmr,vmr,c];

%% Hill 


[obj,grad,Hess] = Hill(Conc,vbs,vEs,vms);


%% Define the symbolic hill function 

syms H(d,b,E,m)

H(d,b,E,m) = log(b + (1-b)/(1+d^m/E));
% Hb = diff(H(d,b,E,m),b);
% HE = diff(H(d,b,E,m),E);
% Hm = diff(H(d,b,E,m),m);
Hgrad = jacobian(H(d,b,E,m),[b,E,m]);
Hhess = jacobian(Hgrad,[b,E,m]);


so = vpa(subs(H(d,b,E,m),[d,b,E,m],[Conc,vbs,vEs,vms]));
sg = vpa(subs(Hgrad,[d,b,E,m],[Conc,vbs,vEs,vms]));
sh = vpa(subs(Hhess,[d,b,E,m],[Conc,vbs,vEs,vms]));

so - obj
sg - grad
sh - Hess

%% Mu function 
tic
[obj,grad,Hess] = mu(Conc,Time(2),[vps,vbetas,vnus,vbs,vEs,vms],'Hess');
t1 = toc

tic

syms mu_syms(t,d,p,beta,nu,b,E,m)
syms mu_syms2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr)

mu_syms(t,d,p,beta,nu,b,E,m) = p*exp(t*(beta-nu+H(d,b,E,m)));
mu_syms2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr) = init*mu_syms(t,d,ps,betas,nus,bs,Es,ms) + init*mu_syms(t,d,pr,betar,nur,br,Er,mr);

% syms mu2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr)
% 
% mu2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr) = mu(t,d,ps,betas,nus,bs,Es,ms) + mu(t,d,pr,betar,nur,br,Er,mr);

mu_grad = jacobian(mu_syms(t,d,p,beta,nu,b,E,m),[p,beta,nu,b,E,m]);
mu_Hess = jacobian(mu_grad,[p,beta,nu,b,E,m]);

so = vpa(subs(mu_syms(t,d,p,beta,nu,b,E,m),[t,d,p,beta,nu,b,E,m],[Time(2),Conc,vps,vbetas,vnus,vbs,vEs,vms]));
sg = vpa(subs(mu_grad,[t,d,p,beta,nu,b,E,m],[Time(2),Conc,vps,vbetas,vnus,vbs,vEs,vms]));
sh = vpa(subs(mu_Hess,[t,d,p,beta,nu,b,E,m],[Time(2),Conc,vps,vbetas,vnus,vbs,vEs,vms]));

t2 = toc

so - obj
sg - squeeze(grad)'
sh - squeeze(Hess)'


%% sig function

tic
[obj,grad,Hess] = sig(Conc,Time(2),[vps,vbetas,vnus,vbs,vEs,vms],'Hess');
t1 = toc

tic
syms lam(d,beta,nu,b,E,m)
syms bn(d,beta,nu,b,E,m)

lam(d,beta,nu,b,E,m) = beta-nu+H(d,b,E,m);
bn(d,beta,nu,b,E,m) = beta+nu-H(d,b,E,m);


syms sig_syms(t,d,p,beta,nu,b,E,m)
syms sig_syms2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr)

sig_syms(t,d,p,beta,nu,b,E,m) = p*bn(d,beta,nu,b,E,m)/lam(d,beta,nu,b,E,m)*(exp(2*lam(d,beta,nu,b,E,m)*t)-exp(t*lam(d,beta,nu,b,E,m)));
sig_syms2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr) ...
= init*sig_syms(t,d,ps,betas,nus,bs,Es,ms) + init*sig_syms(t,d,pr,betar,nur,br,Er,mr);

sig_grad = jacobian(sig_syms(t,d,p,beta,nu,b,E,m),[p,beta,nu,b,E,m]);
sig_Hess = jacobian(sig_grad,[p,beta,nu,b,E,m]);

so = vpa(subs(sig_syms(t,d,p,beta,nu,b,E,m),[t,d,p,beta,nu,b,E,m],[Time(2),Conc,vps,vbetas,vnus,vbs,vEs,vms]));
sg = vpa(subs(sig_grad,[t,d,p,beta,nu,b,E,m],[Time(2),Conc,vps,vbetas,vnus,vbs,vEs,vms]));
sh = vpa(subs(sig_Hess,[t,d,p,beta,nu,b,E,m],[Time(2),Conc,vps,vbetas,vnus,vbs,vEs,vms]));

t2 = toc

so - obj
sg - squeeze(grad)
sh - squeeze(Hess)


%% The objective

tic
[obj,grad,Hess] = BD_Obj_full(theta_c,DATA,Conc,Time,init,'Hess');
t1 = toc

tic
syms obj1(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr)
syms obj2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr)
syms obj_syms(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr)

obj1(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr) ...
= (DATA(2) - mu_syms2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr))^2 ...
             /(2*sig_syms2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr));

obj2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr) ...
= log(1/sqrt(2*pi*sig_syms2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr)));


obj_syms(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr) = obj1(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr) - obj2(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr);

obj_grad = jacobian(obj_syms(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr),[ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr]);
obj_Hess = jacobian(obj_grad,[ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr]);
t2 = toc

tic

so = vpa(subs(obj_syms(t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr),[t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr],[Time(2),Conc,theta_c(1:12)]));
sg = vpa(subs(obj_grad,[t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr],[Time(2),Conc,theta_c(1:12)]));
sh = vpa(subs(obj_Hess,[t,d,ps,betas,nus,bs,Es,ms,pr,betar,nur,br,Er,mr],[Time(2),Conc,theta_c(1:12)]));


t3 = toc


so - obj
sg - squeeze(grad(1:12))
sh - squeeze(Hess(1:12,1:12))


