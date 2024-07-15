%% Load data

load('Result\mix_logistic_S2.mat')

%% Computational time

t = tiledlayout(1,3);
ax = nexttile;

vec_name = {'SQP','SQP grad','IP','IP hess','CRNAS'};
vec_mat  = [fmin_sqp_t_hist,fmin_sqp_grad_t_hist,fmin_IP_t_hist,fmin_IP_hess_t_hist,cubic_t_hist];
plot_vecs(vec_mat,vec_name, [0,1.2*max(max(vec_mat))], 'Computational time', 'Comparison of computational time',[1,0,0;0,1,0;1,1,0;0,0,1;0,1,1],'rank sum',0)

% max(vec_mat)
% ax.YScale = 'linear';


%% Function value
ax = nexttile;
vec_mat  = [fmin_sqp_val_hist,fmin_sqp_grad_val_hist,fmin_IP_val_hist,fmin_IP_hess_val_hist,cubic_val_hist];
plot_vecs(vec_mat,vec_name, [0,1], 'Function value', 'Comparison of function value',[1,0,0;0,1,0;1,1,0;0,0,1;0,1,1],'rank sum',1)
ac = gca;
ac.YLim = [0,10]

%% Iteration
ax = nexttile;
vec_mat  = [fmin_sqp_iter_hist,fmin_sqp_grad_iter_hist,fmin_IP_iter_hist,fmin_IP_hess_iter_hist,cubic_iter_hist];
plot_vecs(vec_mat,vec_name, [0,1.2*max(max(vec_mat))], 'Iteration number', 'Comparison of iteration',[1,0,0;0,1,0;1,1,0;0,0,1;0,1,1],'rank sum',0)


