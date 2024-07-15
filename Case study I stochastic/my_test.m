rng("default")

Time = 0:3:36;
Conc = [0,0.0313,0.0625,0.125,0.25,0.375,0.5,1.25,2.5,3.75,5];  
NT   = length(Time);
ND   = length(Conc);
NR   = 13;
ps = 0.4;
pr = 1-ps;
birth_s = 0.56;
death_s = 0.54;
b_s = 0.87;
E_s = 0.08;
m_s = 2.5;
birth_r = 0.13;
death_r = 0.08;
b_r = 0.97;
E_r = 1;
m_r = 1.5;
c   = 0;
theta = [ps,birth_s,death_s,b_s,E_s,m_s,pr,birth_r,death_r,b_r,E_r,m_r,c];
PMAT  = [ps,birth_s,death_s,b_s,E_s^m_s,m_s;
        pr,birth_r,death_r,b_r,E_r^m_r,m_r]; % Changed of variables

n = 3291;
DATA = sto_gen_bd(NR,Conc,Time,n,theta);

theta_c = [ps,birth_s,death_s,b_s,E_s^m_s,m_s,pr,birth_r,death_r,b_r,E_r^m_r,m_r,c];



%% Grad

my_fun = @(x) BD_Obj_full(x,DATA,Conc,Time,n,'Obj');
my_grad = gradient(my_fun)

%% Time comparison

tic
obj_old = BD_Obj_full(theta_c,DATA,Conc,Time,n,'Obj');
t1 = toc

% tic
% obj_old = BD_Obj_full(theta_c,DATA,Conc,Time,n,'grad');
% t2 = toc

tic
[obj,grad,Hess] = BD_Obj_full(theta_c,DATA,Conc,Time,n,'Hess');
t3 = toc

tic
[obj_f,grad_f,Hess_f] = BD_Obj_full_finite(theta_c,DATA,Conc,Time,n,1e-10,1e-6);
t4 = toc


%%  Check for the Obj gradient computation


eps = 1e-11;


v1   = zeros(1,13);
i = 5;
v1(i) = 1;
v2  = zeros(1,13);
j = 6;
v2(j) = 1;





tic
[obj1,grad1,Hess1] = BD_Obj_full(theta_c,DATA,Conc,Time,n,'Hess');
t1 = toc
tic
[obj2,grad2,Hess2] = BD_Obj_full_finite(theta_c,DATA,Conc,Time,n,1e-10,1e-7);
t2 = toc
% [obj2,~,~] = BD_Obj_full(theta_c+eps*v1,DATA,Conc,Time,n,'Hess');
% [obj3,~,~] = BD_Obj_full(theta_c+eps*v2,DATA,Conc,Time,n,'Hess');
% [obj4,~,~] = BD_Obj_full(theta_c+eps*(v1+v2),DATA,Conc,Time,n,'Hess');

% diff_grad = ((obj2 - obj)./eps - grad(i))./grad(i)
% diff_Hess = ((obj4 - obj2 - obj3 + obj)./(eps^2) - Hess(i,j))./Hess(i,j)

diff_grad = norm(grad1 - grad2)./norm(grad1)
diff_Hess = norm(Hess1-Hess2)./norm(Hess1)


%% Random check

eps = 1e-7;

d = rand(1,13);
d = d./norm(d);

[obj2,~,~] = BD_Obj_full(theta_c+eps*d,DATA,Conc,Time,n,'grad');

[obj3,~,~] = BD_Obj_full(theta_c-eps*d,DATA,Conc,Time,n,'grad');

grad_diff = ((obj2 - obj)./eps - d*grad')./(d*grad')
Hess_diff = ((obj2 + obj3 - 2*obj )./(eps^2) - d*Hess*d')./(d*Hess*d')



%% Check for the objective computation

DATA_old = permute(DATA,[3,2,1]);

tic
obj_old = BD_Obj(theta,DATA_old,Conc,Time,NR,2)
t1 = toc
tic
[obj_new,~,~]= BD_Obj_full(theta_c,DATA,Conc,Time,n)
t2 = toc


%% Check for the Hill gradient computation


eps = 1e-6;

vec = [b_s,E_s^m_s,m_s];
d = rand(1,3);
d = d./norm(d);
vec_eps = vec + eps*d;
b = vec(1);
E = vec(2);
m = vec(3);

c_num = 3;

[obj_H,grad_H,Hess_H] = Hill(Conc(c_num), b,E,m);

[obj2,~,~] = Hill(Conc(c_num),b_s + eps*d(1),E_s^m_s + eps*d(2), m_s + eps*d(3));
[obj3,~,~] = Hill(Conc(c_num),b_s - eps*d(1),E_s^m_s - eps*d(2), m_s - eps*d(3));

grad_diff = (obj2 - obj_H)./eps - d*grad_H'
Hess_diff = (obj2 + obj3 - 2*obj_H )./(eps^2) - d*Hess_H*d'



% [obj_bl,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s);
% [obj_br,~,~] = Hill(Conc(c_num), b_s+eps,E_s^m_s,m_s);
% [obj_El,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s);
% [obj_Er,~,~] = Hill(Conc(c_num), b_s,E_s^m_s+eps,m_s);
% [obj_ml,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s);
% [obj_mr,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s+eps);
% 
% [(obj_br-obj_bl)/(eps) - grad_H(1),(obj_Er-obj_El)/(eps) - grad_H(2),(obj_mr-obj_ml)/(eps) - grad_H(3)]



% [obj_bb,~,~]  = Hill(Conc(c_num), b_s+2*eps,E_s^m_s,m_s);
% [obj_EE,~,~]  = Hill(Conc(c_num), b_s,E_s^m_s+2*eps,m_s);
% [obj_mm,~,~]  = Hill(Conc(c_num), b_s,E_s^m_s,m_s+2*eps);
% [obj_bE1,~,~] = Hill(Conc(c_num), b_s+eps,E_s^m_s+eps,m_s);
% [obj_bE2,~,~] = Hill(Conc(c_num), b_s+eps,E_s^m_s,m_s);
% [obj_bE3,~,~] = Hill(Conc(c_num), b_s,E_s^m_s+eps,m_s);
% [obj_bE4,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s);
% [obj_bm1,~,~] = Hill(Conc(c_num), b_s+eps,E_s^m_s,m_s+eps);
% [obj_bm2,~,~] = Hill(Conc(c_num), b_s+eps,E_s^m_s,m_s);
% [obj_bm3,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s+eps);
% [obj_bm4,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s);
% [obj_Em1,~,~] = Hill(Conc(c_num), b_s,E_s^m_s+eps,m_s);
% [obj_Em2,~,~] = Hill(Conc(c_num), b_s,E_s^m_s+eps,m_s+eps);
% [obj_Em3,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s);
% [obj_Em4,~,~] = Hill(Conc(c_num), b_s,E_s^m_s,m_s+eps);
% 
% Hess_diff = zeros(3,3);
% Hess_diff(1,1) = (obj_bb - 2*obj_bE2+obj_bE4)/(eps^2);
% Hess_diff(2,2) = (obj_EE - 2*obj_Em1+obj_Em3)/(eps^2);
% Hess_diff(3,3) = (obj_mm - 2*obj_Em4+obj_Em3)/(eps^2);
% Hess_diff(1,2) = (obj_bE1 - obj_bE2 - obj_bE3 + obj_bE4)/(eps^2);
% Hess_diff(1,3) = (obj_bm1 - obj_bm2 - obj_bm3 + obj_bm4)/(eps^2);
% Hess_diff(2,3) = (obj_Em2 - obj_Em1 - obj_Em4 + obj_Em3)/(eps^2);
% Hess_diff=Hess_diff+Hess_diff.'-diag(diag(Hess_diff));
% 
% Hess_diff

%% Check for the mu gradient computation

eps = 1e-8;
d = rand(1,6);
d = d./norm(d);


p    = pr;
beta = birth_r;
nu   = death_r;
b    = b_r;
E    = E_r;
m    = m_r;
vec = [p,beta,nu,b,E,m];

i = 2;
j = 4;

[obj,grad,Hess] = mu(Conc,Time',vec,'Hess');
[obj2,~,~] = mu(Conc,Time',vec+eps*d,'grad');
[obj3,~,~] = mu(Conc,Time',vec-eps*d,'grad');

grad_diff = (obj2(i,j) - obj(i,j))./eps - d*squeeze(grad(i,j,:))
Hess_diff = (obj2(i,j) + obj3(i,j) - 2*obj(i,j) )./(eps^2) - d*squeeze(Hess(i,j,:,:))*d'


%% Check for the ratio gradient computation
eps = 1e-5;
d = rand(1,5);
d = d./norm(d);


p    = pr;
beta = birth_r;
nu   = death_r;
b    = b_r;
E    = E_r;
m    = m_r;


i = 7;
j = 4;

[obj,grad,Hess] = ratio(Conc(i),beta,nu,b,E,m);
[obj2,~,~] = ratio(Conc(i),beta+eps*d(1),nu+eps*d(2),b+eps*d(3),E+eps*d(4),m+eps*d(5));
[obj3,~,~] = ratio(Conc(i),beta-eps*d(1),nu-eps*d(2),b-eps*d(3),E-eps*d(4),m-eps*d(5));

grad_diff = (obj2 - obj)./eps - d*grad'
Hess_diff = (obj2 + obj3 - 2*obj )./(eps^2) - d*Hess*d'





%%  Check for the sig gradient computation ************
eps = 1e-5;
d = rand(1,6);
d = d./norm(d);


p    = pr;
beta = birth_r;
nu   = death_r;
b    = b_r;
E    = E_r;
m    = m_r;
vec = [p,beta,nu,b,E,m];


i = 2;
j = 2;

[obj,grad,Hess] = sig(Conc(i),Time(j),vec,'Hess');
[obj2,~,~] = sig(Conc(i),Time(j),vec+eps*d,'grad');
[obj3,~,~] = sig(Conc(i),Time(j),vec-eps*d,'grad');

grad_diff = (obj2 - obj)./eps - d*grad'
Hess_diff = (obj2 + obj3 - 2*obj )./(eps^2) - d*Hess*d'








%%
eps = 1e-5;

obj_a = sig(Conc(i),Time(j),p+eps,beta,nu,b,E,m);
obj_b = sig(Conc(i),Time(j),p,beta,nu,b,E+eps,m);
obj_ab = sig(Conc(i),Time(j),p+eps,beta,nu,b,E+eps,m);

diff = (obj_ab - obj_a -obj_b+obj)./(eps^2) - Hess(1,5)
