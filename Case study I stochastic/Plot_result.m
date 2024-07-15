%% load data

load('Result\BD_opt20_100_uc.mat')

%% Computational time
t = tiledlayout(1,3);
ax = nexttile;

vec_name = {'IP hess','SQP grad','CRNAS'};
vec_mat  = [t_IP_hist',t_sqp_hess_hist',t_cubic_hist'];
plot_vecs(vec_mat,vec_name, [0,max(max(vec_mat))], 'Computational time', 'Comparison of computational time',[1,0,0;0,1,0;0,0,1],'rank sum',0)

ax = gca;
% ax.YScale = 'log';
ax.YLim = [0,1.3*max(max(vec_mat))];

%% Function value

ax = nexttile;

vec_mat  = [modified_obj_IP_hist',modified_obj_sqp_hess_hist',modified_obj_cubic_hist'];
% vec_mat(31:32,:) = []
plot_vecs(vec_mat,vec_name, [0,max(max(vec_mat))], 'Objective value', 'Comparison of objective value',[1,0,0;0,1,0;0,0,1],'rank sum',1)
ac = gca;
% ac.YScale = 'log';
ac.YLim   = [min(min(vec_mat)),1.1*max(max(vec_mat))];


%% Iteration 


ax = nexttile;

vec_mat  = [iter_IP_hist',iter_sqp_hess_hist',iter_cubic_hist'];
% vec_mat(31:32,:) = []
plot_vecs(vec_mat,vec_name, [0,max(max(vec_mat))], 'Iteration number', 'Comparison of iteration',[1,0,0;0,1,0;0,0,1],'rank sum',0)
ac = gca;
% ac.YScale = 'log';
ac.YLim = [0,1.3*max(max(vec_mat))];

