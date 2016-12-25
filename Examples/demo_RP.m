
clear
startup
randn('state',0); rand('state',0);


ndim  = [5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3];
rrank = [50,  50,  50, 50,  50,  50, 50,  50,  50]; 
pfac  = [1.5, 2.0, 2.3, 2.5, 2.8, 3.0, 3.3, 3.5, 3.8];

%
%noise strength; may change from 0.0001 to 0.05
sigma = 0.001; 

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
S = diag(1000*(abs(randn(rankk,1)).^2)); %kai_square
L = L0*S;
dof = rankk*(m+n-rankk);
rank_reconstruct = rankk;
X = L*R';

% make random sampling, problem and initial guess
samples = floor(OS * dof);
Omega = make_rand_Omega(m,n,samples);
prob = make_prob(L,R,Omega,rank_reconstruct,sigma); % <- you can choose another rank here

%--------------------------------------------------------------------------
%LRGeomCG with rank = 1.00*rankk,
%--------------------------------------------------------------------------
tstart = clock;
rank_est = round(1.00*rankk);
prob.r = rank_est;
options = default_opts();
options.maxit = 150;
x0 = make_start_x(prob);
prob.t0 = tstart;
prob.norm_M_Omega = norm(prob.data);
[Xcg_CG,hist_CG] = LRGeomCG(prob,options,x0);
out_time_LRGeomCG = etime(clock,tstart);
Xcg2 = Xcg_CG.U * diag(Xcg_CG.sigma) * Xcg_CG.V';
error_cg = norm(Xcg2-X,'fro')/norm(X,'fro')

%--------------------------------------------------------------------------
%RP
%--------------------------------------------------------------------------

opts = default_opts_RP();

%since this example is a noiseless case, we use small opts.rel_out_obj_tol
%and opts.rel_out_obj_diff_tol
% tollerance w.r.t. objective values
opts.rel_out_obj_tol = 0.0001; 
% tollerance w.r.t. objective value difference
opts.rel_out_obj_diff_tol = 0.000001; 


tstart = clock;
[Xcg_RP, histout_RP, hist_out_all] = LRGeomCG_RP(prob, opts);
out_time_RP = etime(clock,tstart);
Xcg2_RP = Xcg_RP.U * diag(Xcg_RP.sigma) * Xcg_RP.V';
error_RP = norm(Xcg2_RP-X,'fro')/norm(X,'fro');

%--------------------------------------------------------------------------
%LRGeomCG with rank detected by RP. 
%--------------------------------------------------------------------------
prob.r = rank(Xcg_RP.U);
options.maxit = 150;
tstart = clock;
prob.t0 = tstart;
x0 = make_start_x(prob);
prob.norm_M_Omega = norm(prob.data);
[Xcg_cg_obj,hist_CG_obj] = LRGeomCG(prob,options,x0);
out_time_CG_obj = etime(clock,tstart)
Xcg2_cg_obj = Xcg_cg_obj.U * diag(Xcg_cg_obj.sigma) * Xcg_cg_obj.V';
error_cg_obj = norm(Xcg2_cg_obj-X,'fro')/norm(X,'fro')


%objective values w.r.t. time slot

figure
semilogy(histout_RP(:,5), histout_RP(:,2)/histout_RP(1,2),'-r','LineWidth',3)
hold on
semilogy(hist_CG(:,5), hist_CG(:,2)/hist_CG(1,2),'-m','LineWidth',3)
hold on
semilogy(hist_CG_obj(:,5), hist_CG_obj(:,2)/hist_CG_obj(1,2),'-','LineWidth',3)

legend('RP','LRGeomCG (1.00r)','LRGeomCG (rank by RP)', 'Location','NE')
axis([0 50 0 1])
grid
box on
set(gca,'FontSize',16); 
xlabel('Time (seconds)','FontSize',16)
ylabel('Relative objective value','FontSize',16);    


%objective values w.r.t. iterations
figure
semilogy(histout_RP(:,2)/histout_RP(1,2),'-r','LineWidth',3)
hold on
semilogy(hist_CG(:,2)/hist_CG(1,2),'-m','LineWidth',3)
hold on
semilogy(hist_CG_obj(:,2)/hist_CG_obj(1,2),'-','LineWidth',3)

legend('RP','LRGeomCG (1.00r)','LRGeomCG (rank by RP)', 'Location','NE')
axis([0 150 0 1])
grid
box on
set(gca,'FontSize',16); 
xlabel('Iterations','FontSize',16)
ylabel('Relative objective value','FontSize',16);  
end




