% =========================================================================
% This script runs the parameter sensitivity expr in IJCAI 2015 paper
%
% Yan, Y., Tan, M., Tsang, I., Yang, Y., Zhang, C., & Shi, Q. (2015).
% Scalable Maximum Margin Matrix Factorization by Active Riemannian Subspace Search. 
% In IJCAI Proceedings-International Joint Conference on Artificial Intelligence.
%
% Contact: yanyan.tju@gmail.com
% =========================================================================

clear
startup
randn('state',0); rand('state',0);

ndim  = [1e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3];
rrank = [20,  50,  50, 50,  50,  50, 50,  50,  50]; 
pfac  = [3.5, 2.0, 2.3, 2.5, 2.8, 3.0, 3.3, 3.5, 3.8];

%
%noise strength; may change from 0.0001 to 0.05
sigma = 0.000; 

for k = 1:1
    
rankk = rrank(k);
OS = pfac(k);
% dimensions
m = ndim(k); n = ndim(k); 
% random factors
L0 = randn(m, rankk); 
[L0,~] = qr(L0,0);
R = randn(n, rankk); 
[R,~] = qr(R,0);
% S = diag(1000*(abs(randn(rankk,1)))); %
S = diag(1000*(abs(randn(rankk,1)))); %kai_square.^2
L = L0*S;
dof = rankk*(m+n-rankk);
rank_reconstruct = rankk;
X = L*R';

X_rating = X;
X_rating(X_rating>0) = 1;
X_rating(X_rating<=0) = -1;

rating_levels = sort(unique(X_rating));

% make random sampling, problem and initial guess
samples = floor(OS * dof);
Omega = make_rand_Omega(m,n,samples);
% prob = make_prob(L,R,Omega,rank_reconstruct,sigma); % <- you can choose another rank here
prob = make_binary_prob(L,R,Omega,rank_reconstruct,sigma); % <- you can choose another rank here

theta = zeros(size(rating_levels,1)-1,1);
for r = 1:size(theta)
    theta(r) = 1/2 * (rating_levels(r) + rating_levels(r+1));
end
prob.theta = theta;
prob.rating_levels = rating_levels;

prob.loss_flag = '1MM';  % binary case: hinge loss
prob.theta_flag = 'fix';  % fix theta
prob.rmse_flag = 1;  % compute the rmse in each iteration

rank_est = round(1.00*rankk);
prob.r = 20;
options = default_opts();
options.maxit = 100;
options.rel_tol_change_res = 1e-8;
x0 = make_start_x(prob);

%%new options for test_rmse
prob.X_test_set = X_rating-prob.temp_omega;
prob.n_test = nnz(prob.X_test_set);
prob.n_train = samples;
prob.train_indices = find(abs(prob.temp_omega)>0);
prob.test_indices = find(abs(prob.X_test_set)>0);
prob.norm_M_Omega = norm(prob.data);

%% mu = 0.2
prob.mu = 0.02;
fprintf('mu=%0.4f\t',prob.mu)
tstart = clock;
prob.t0 = tstart;
[Xcg_CG,hist_CG_02] = LRGeomCG_bin(prob,options,x0);
fprintf('done...\n')

%% mu = 0.01
prob.mu = 0.01;
fprintf('mu=%0.4f\t',prob.mu)
tstart = clock;
prob.t0 = tstart;
[Xcg_CG,hist_CG_01] = LRGeomCG_bin(prob,options,x0);
fprintf('done...\n')

%% mu = 0.005
prob.mu = 0.005;
fprintf('mu=%0.4f\t',prob.mu)
tstart = clock;
prob.t0 = tstart;
[Xcg_CG,hist_CG_005] = LRGeomCG_bin(prob,options,x0);
fprintf('done...\n')

%% mu = 0.001
prob.mu = 0.001;
fprintf('mu=%0.4f\t',prob.mu)
tstart = clock;
prob.t0 = tstart;
[Xcg_CG,hist_CG_001] = LRGeomCG_bin(prob,options,x0);
fprintf('done...\n')

%% mu = 0.0001
prob.mu = 0.0001;
fprintf('mu=%0.4f\t',prob.mu)
tstart = clock;
prob.t0 = tstart;
[Xcg_CG,hist_CG_0001] = LRGeomCG_bin(prob,options,x0);
fprintf('done...\n')

%% mu = 0.0000
prob.mu = 0.0000;
fprintf('mu=%0.4f\t',prob.mu)
tstart = clock;
prob.t0 = tstart;
[Xcg_CG,hist_CG_0] = LRGeomCG_bin(prob,options,x0);
fprintf('done...\n')

%%
lineWidth = 3;
%% RMSE on testing points
font_size = 21;
figure
plot([1:length(hist_CG_02(:,8))], hist_CG_02(:,8),'-k','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_01(:,8))], hist_CG_01(:,8),'-r','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_005(:,8))], hist_CG_005(:,8),'-c','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_001(:,8))], hist_CG_001(:,8),'-.g','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_0001(:,8))], hist_CG_0001(:,8),'-.m','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_0(:,8))], hist_CG_0(:,8),'-.b','LineWidth',lineWidth); hold on

grid
box on
legend('\lambda = 0.02','\lambda = 0.01','\lambda = 0.005','\lambda = 0.001','\lambda = 0.0001','\lambda = 0')
axis([0 100 0.755 0.8])
set(gca,'FontSize',font_size-5); 
xlabel('# iterations','FontSize',font_size)
ylabel('Testing RMSE','FontSize',font_size);    

%% RMSE on testing points (another method to compute RMSE. The result is identical with the above result)

figure

plot([1:length(hist_CG_02(:,9))], hist_CG_02(:,9),'-k','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_01(:,9))], hist_CG_01(:,9),'-r','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_005(:,9))], hist_CG_005(:,9),'-c','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_001(:,9))], hist_CG_001(:,9),'-.g','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_0001(:,9))], hist_CG_0001(:,9),'-.m','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_0(:,9))], hist_CG_0(:,9),'-.b','LineWidth',lineWidth); hold on

grid
box on
legend('\lambda = 0.02','\lambda = 0.01','\lambda = 0.005','\lambda = 0.001','\lambda = 0.0001','\lambda = 0')
axis([0 100 0.755 0.8])
set(gca,'FontSize',font_size-5); 
xlabel('# iterations','FontSize',font_size)
ylabel('Testing RMSE','FontSize',font_size);     

%% RMSE on training points
figure

plot([1:length(hist_CG_02(:,7))], hist_CG_02(:,7),'-k','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_01(:,7))], hist_CG_01(:,7),'-r','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_005(:,7))], hist_CG_005(:,7),'-c','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_001(:,7))], hist_CG_001(:,7),'-.g','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_0001(:,7))], hist_CG_0001(:,7),'-.m','LineWidth',lineWidth); hold on

plot([1:length(hist_CG_0(:,7))], hist_CG_0(:,7),'-.b','LineWidth',lineWidth); hold on

grid
box on
legend('\lambda = 0.02','\lambda = 0.01','\lambda = 0.005','\lambda = 0.001','\lambda = 0.0001','\lambda = 0')
set(gca,'FontSize',font_size-5); 
xlabel('# iterations','FontSize',font_size)
ylabel('Training RMSE','FontSize',font_size);     


%% objetive value 
figure
semilogy([1:length(hist_CG_02(:,7))], hist_CG_02(:,2)/hist_CG_02(1,2),'-k','LineWidth',lineWidth); hold on

hold on
semilogy([1:length(hist_CG_01(:,7))], hist_CG_01(:,2)/hist_CG_01(1,2),'-r','LineWidth',lineWidth)

hold on
semilogy([1:length(hist_CG_005(:,7))], hist_CG_005(:,2)/hist_CG_005(1,2),'-c','LineWidth',lineWidth)

hold on
semilogy([1:length(hist_CG_001(:,7))], hist_CG_001(:,2)/hist_CG_001(1,2),'-g','LineWidth',lineWidth)

hold on
semilogy([1:length(hist_CG_0001(:,7))], hist_CG_0001(:,2)/hist_CG_0001(1,2),'-m','LineWidth',lineWidth)

hold on
semilogy([1:length(hist_CG_0(:,7))], hist_CG_0(:,2)/hist_CG_0(1,2),'-b','LineWidth',lineWidth)

grid
box on
legend('\lambda = 0.02','\lambda = 0.01','\lambda = 0.005','\lambda = 0.001','\lambda = 0.0001','\lambda = 0')
axis([0 100 0 1])
set(gca,'FontSize',font_size-5); 
xlabel('# iterations','FontSize',font_size)
ylabel('Relative Objective Value','FontSize',font_size);


end




