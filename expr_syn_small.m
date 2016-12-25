% =========================================================================
% This script runs the small synthetic expr in IJCAI 2015 paper
%
% Yan, Y., Tan, M., Tsang, I., Yang, Y., Zhang, C., & Shi, Q. (2015).
% Scalable Maximum Margin Matrix Factorization by Active Riemannian Subspace Search. 
% In IJCAI Proceedings-International Joint Conference on Artificial Intelligence.
%
% Contact: yanyan.tju@gmail.com
% =========================================================================

%%
% clear
%%
% if isempty(strcmp(computer, 'WIN'))
% matlabpool('open',4)
% end
%% add path
addpath(genpath('./'));

%% set the choices
do_RP_M3F = 1;
do_RP_LS = 1;
do_CG_M3F = 1;
do_CG_LS = 1;

prefix = 'expr_data_';
expr_No = input('Choose a synthetic experiment data No.:', 's');

result_root = ['small_syn_result', expr_No];
if ~exist(result_root, 'dir')
    system(['mkdir ', result_root]);
end

%%
extra_est_rank = 0;

rmse_flag = 0;
%% generate synthetic data

% set the directory to store the generated data
data_save_path = [prefix, expr_No];
if exist([prefix, expr_No])
    display('This direxctory exists. Loading...')
    data_save_file = [data_save_path, '/synthetic_data'];
    load(data_save_file)
    threshold = prob.theta;
else
    mkdir([prefix, expr_No])
    mkdir([prefix, expr_No, '/result'])

    m = 1000;  % # of rows
    n = 1000;  % # of cols
    data_rank = 20;  % rank
    OS = 3.5;

    rating_levels = [1;2;3;4;5];
    theta = zeros(size(rating_levels,1)-1,1);
    for r = 1:size(theta)
        theta(r) = 1/2 * (rating_levels(r) + rating_levels(r+1));
    end

    %noise strength; may change from 0.0001 to 0.05
    sigma = 0.000; 

    % startup
    scale = 1;
    bias = 0;
    L0 = randn(m, data_rank); 
    [L0,~] = qr(L0,0);
    R = randn(n, data_rank); 
    [R,~] = qr(R,0);
    % S = diag(1000*(abs(randn(rank,1)))); %
    S = diag(1000*(abs(randn(data_rank,1)))); %kai_square
    L = scale*L0*S;
%     L = L0*S;
    % =========================================================================
    dof = data_rank*(m+n-data_rank);
    rank_reconstruct = data_rank;
    X =  L*R'+bias;

        % make random sampling, problem and initial guess

    samples = floor(OS * dof);
    Omega = make_rand_Omega(m, n, samples);

%     prob = make_prob_rating(L, R, Omega, rank_reconstruct, sigma, rating_levels, theta, bias); % <- you can choose another rank here
    prob = make_rating_prob_bias(L, R, Omega, rank_reconstruct, sigma, rating_levels, bias); % <- you can choose another rank here
    prob.theta = theta;
    prob.rating_levels = rating_levels;

    % Reconstruct the observed matrix (y)
    X_observed = sparse(prob.Omega_i, prob.Omega_j, prob.data, prob.n1, prob.n2 );

    % Output the observed matrix as tsv files
    [i, j, val] = find(X_observed);
    data_dump = [i-1, j-1, val];
    fout = fopen([data_save_path, '/prob_rating_hog.tsv'],'w');
    fprintf(fout, '%d\t%d\t%d\n', data_dump');
    fclose(fout);

    % Reconstruct the original matrix X
        % discrete
        cuts = prob.cuts;
        X_rating = ones(size(X));
        num_levels = size(rating_levels,1);
        for i=1:num_levels-1
            X_rating(X>cuts(i)) = i+1;
        end

    % Output the original matrix (ground truth) as tsv file
    [i, j, val] = find(X_rating);  % make sure that there is no '0' in the data
    data_dump = [i-1, j-1, val];
    fout = fopen([data_save_path, '/org_rating_hog.tsv'],'w');
    fprintf(fout, '%d\t%d\t%d\n', data_dump');
    fclose(fout);

    % Output the testing matrix as tsv file
    X_test = X_rating - X_observed;
    [i, j, val] = find(X_test);
    data_dump = [i-1, j-1, val];
    fout = fopen([data_save_path, '/test_rating_hog.tsv'],'w');
    fprintf(fout, '%d\t%d\t%d\n', data_dump');
    fclose(fout);

    data_save_file = [data_save_path, '/synthetic_data'];
    save(data_save_file, 'X_rating', 'X_test', 'prob', 'X_observed');

end

prob.train_indices = prob.Omega;
prob.test_indices = setdiff(1:size(X_rating,1)*size(X_rating,2), prob.train_indices);
Known_train = prob.train_indices;
Known_test = prob.test_indices;
prob.X_test_mat = X_rating-prob.temp_omega;    

CG_est_rank = prob.r;  % the estimation of matrix rank for CG-LS and CG-MM
GROUSE_est_rank = prob.r;
SCGMC_est_rank = prob.r;

%%
% ====================================================
% ------------------------ RP ------------------------
% ====================================================
%% RP-M3F
if do_RP_M3F

    display('START RP-M3F')
    prob.loss_flag = 'RMM';
     prob.theta_flag = 'update';

    opts = default_opts_RP();
    opts.tar_maxit = 60;

    %since this example is a noiseless case, we use small opts.rel_out_obj_tol
    %and opts.rel_out_obj_diff_tol
    % tollerance w.r.t. objective values
    opts.rel_out_obj_tol = 0.175;
    opts.rel_out_obj_diff_tol = 0.000001; 
    opts.eta = 0.75;
    opts.tar_rel_inner_tol = 0.0001;
    opts.rel_inner_tol = 0.5;
    
    theta = prob.theta;

    prob.mu = 0.000;
    prob.rmse_flag = rmse_flag;
    prob.max_iter_out = 50;

    start_time = clock;
    prob.t0 = start_time;

    [Xcg_RPMM_theta, hist_RPMM_theta, theta_RPMM_theta, hist_out_all_RPMM_theta] = LRGeomCG_RP_mix(prob, opts, theta);

    out_time_RPMM = etime(clock, start_time);
    Xcg2_RPMM_theta = Xcg_RPMM_theta.U * diag(Xcg_RPMM_theta.sigma) * Xcg_RPMM_theta.V';
    rank_RPMM_theta = rank(Xcg_RPMM_theta.U);

    % discrete data
    cuts = theta_RPMM_theta;
    X_RPMM_predict_theta = ones(size(Xcg2_RPMM_theta));
    rating_levels = prob.rating_levels;
    for i=1:size(rating_levels,1)-1
        X_RPMM_predict_theta(Xcg2_RPMM_theta>cuts(i)) = i+1;
    end

    rte_test_RPMM = norm(X_RPMM_predict_theta(Known_test)-X_test(Known_test),'fro')/norm(X_test(Known_test),'fro');
    rte_train_RPMM = norm(X_RPMM_predict_theta(Known_train)-X_observed(Known_train),'fro')/norm(X_observed(Known_train),'fro');

    rmse_test_RPMM_rating = norm(X_RPMM_predict_theta(Known_test) - X_test(Known_test), 'fro') / sqrt(length(Known_test));
    rmse_train_RPMM_rating = norm(X_RPMM_predict_theta(Known_train) - X_observed(Known_train), 'fro') / sqrt(length(Known_train));

    rmse_test_RPMM_org = norm(Xcg2_RPMM_theta(Known_test) - X_test(Known_test), 'fro') / sqrt(length(Known_test));
    rmse_train_RPMM_org = norm(Xcg2_RPMM_theta(Known_train) - X_observed(Known_train), 'fro') / sqrt(length(Known_train));

    save_file = [result_root, '/RP_RMM'];
    save(save_file, 'hist_RPMM_theta', 'out_time_RPMM', 'rte_test_RPMM', 'rte_train_RPMM', 'rmse_test_RPMM_rating', 'rmse_train_RPMM_rating', 'rmse_test_RPMM_org', 'rmse_train_RPMM_org', 'theta_RPMM_theta', 'rank_RPMM_theta')

end

%% RP_LS
if do_RP_LS

    display('START RP-LS')
    prob.loss_flag = 'LS';
    prob.theta_flag = 'fix';

    opts = default_opts_RP();

    opts.rel_out_obj_tol = 0.09; 
    % tollerance w.r.t. objective value difference
    opts.rel_out_obj_diff_tol = 0.00001; 
    opts.eta = 0.65;
    opts.tar_maxit = 80;
    opts.tar_rel_inner_tol = 0.0000001;

    prob.rmse_flag = rmse_flag;
    prob.mu = 0.000;
    prob.max_iter_out = 50;

    start_time = clock;
    theta = prob.theta;
    [Xcg_RPLS_theta, hist_RPLS_theta, theta_RPLS_theta, hist_out_all_RPLS_theta] = LRGeomCG_RP_mix(prob, opts, theta);

    out_time_RPLS = etime(clock, start_time);
    Xcg2_RPLS_theta = Xcg_RPLS_theta.U * diag(Xcg_RPLS_theta.sigma) * Xcg_RPLS_theta.V';
    rank_RPLS_theta = rank(Xcg_RPLS_theta.U);

    % discrete data
    cuts = theta_RPLS_theta;
    X_RPLS_predict_theta = ones(size(Xcg2_RPLS_theta));
    rating_levels = prob.rating_levels;
    for i=1:size(rating_levels,1)-1
        X_RPLS_predict_theta(Xcg2_RPLS_theta>cuts(i)) = i+1;
    end

    rmse_test_RPLS_org = sqrt(sum((Xcg2_RPLS_theta(Known_test) - X_test(Known_test)).^2) / length(Known_test));  
    rmse_train_RPLS_org = sqrt(sum((Xcg2_RPLS_theta(Known_train) - X_observed(Known_train)).^2) / length(Known_train));

    rmse_test_RPLS_rating = sqrt(sum((X_RPLS_predict_theta(Known_test) - X_test(Known_test)).^2) / length(Known_test));  
    rmse_train_RPLS_rating = sqrt(sum((X_RPLS_predict_theta(Known_train) - X_observed(Known_train)).^2) / length(Known_train));

    rte_test_RPLS = norm(X_RPLS_predict_theta(Known_test)-X_test(Known_test),'fro')/norm(X_test(Known_test),'fro');
    rte_train_RPLS = norm(X_RPLS_predict_theta(Known_train)-X_observed(Known_train),'fro')/norm(X_observed(Known_train),'fro');

    save_file = [result_root, '/RP_LS'];
    save(save_file, 'hist_RPLS_theta', 'out_time_RPLS', 'theta_RPLS_theta', 'rte_test_RPLS', 'rte_train_RPLS', 'rank_RPLS_theta', 'rmse_test_RPLS_rating', 'rmse_train_RPLS_rating', 'rmse_test_RPLS_org', 'rmse_train_RPLS_org')

end
%% CG_M3F
if do_CG_M3F
%     ss = input('', s);
%     rank_est_CGMM = CG_est_rank;
    rank_est_CGMM = prob.r;
    display('START CG-M3F')

    prob.loss_flag = 'RMM';
    prob.theta_flag = 'update';

    options = default_opts();
    options.maxit = 100;
    x0 = make_start_x(prob);
    prob.norm_M_Omega = norm(prob.data);

    rank_bk = prob.r;
    prob.r = rank_est_CGMM + extra_est_rank;
    prob.rmse_flag = rmse_flag;
    prob.mu = 0.001;

    start_time = clock;
    prob.t0 = start_time;

    [Xcg_CGMM_theta, hist_CGMM_theta, theta_CGMM_theta] = LRGeomCG_mix(prob,options,x0);

    out_time_CGMM = etime(clock,start_time);

    prob.r = rank_bk;

    Xcg2_theta = Xcg_CGMM_theta.U * diag(Xcg_CGMM_theta.sigma) * Xcg_CGMM_theta.V';

    % discrete data
    cuts = theta_CGMM_theta;
    X_CGMM_predict_theta = ones(size(Xcg2_theta));
    rating_levels = prob.rating_levels;
    for i=1:size(rating_levels,1)-1
        X_CGMM_predict_theta(Xcg2_theta>cuts(i)) = i+1;
    end

    rte_test_CGMM = norm(X_CGMM_predict_theta(Known_test)-X_test(Known_test),'fro')/norm(X_test(Known_test),'fro');
    rte_train_CGMM = norm(X_CGMM_predict_theta(Known_train)-X_observed(Known_train),'fro')/norm(X_observed(Known_train),'fro');

    rmse_test_CGMM_rating = norm(X_CGMM_predict_theta(Known_test) - X_test(Known_test), 'fro') / sqrt(length(Known_test));
    rmse_train_CGMM_rating = norm(X_CGMM_predict_theta(Known_train) - X_observed(Known_train), 'fro') / sqrt(length(Known_train));

    rmse_test_CGMM_org = norm(Xcg2_theta(Known_test) - X_test(Known_test), 'fro') / sqrt(length(Known_test));
    rmse_train_CGMM_org = norm(Xcg2_theta(Known_train) - X_observed(Known_train), 'fro') / sqrt(length(Known_train));

    save_file = [result_root, '/CG_RMM'];
    save(save_file,'hist_CGMM_theta', 'out_time_CGMM', 'rte_test_CGMM', 'rte_train_CGMM', 'rmse_test_CGMM_rating', 'rmse_train_CGMM_rating', 'rmse_test_CGMM_org', 'rmse_train_CGMM_org')
end

%% CG_LS

if do_CG_LS
    rank_est_CGLS = rank_RPMM_theta;
    display('START CG-LS')

    prob.loss_flag = 'LS';
    prob.theta_flag = 'fix';

    options = default_opts();
    options.maxit = 300;
    x0 = make_start_x(prob);
    prob.norm_M_Omega = norm(prob.data);

    rank_bk = prob.r;
    prob.r = rank_est_CGLS + extra_est_rank;
    prob.rmse_flag = rmse_flag;

    start_time = clock;
    prob.t0 = start_time;

    [Xcg_CGLS_theta, hist_CGLS_theta, theta_CGLS_theta] = LRGeomCG_mix(prob,options,x0);  % input variables: LRGeomCG(prob, opts, xc_new_init, fc_new_init, rel_inner_tol, theta, flag);
    out_time_CGLS = etime(clock,start_time);

    prob.r = rank_bk;

    Xcg2_theta = Xcg_CGLS_theta.U * diag(Xcg_CGLS_theta.sigma) * Xcg_CGLS_theta.V';

    % discrete data
    cuts = theta_CGLS_theta;
    X_CGLS_predict_theta = ones(size(Xcg2_theta));
    rating_levels = prob.rating_levels;
    for i=1:size(rating_levels,1)-1
        X_CGLS_predict_theta(Xcg2_theta>cuts(i)) = i+1;
    end

    rte_test_CGLS = norm(X_CGLS_predict_theta(Known_test)-X_test(Known_test),'fro')/norm(X_test(Known_test),'fro');
    rte_train_CGLS = norm(X_CGLS_predict_theta(Known_train)-X_observed(Known_train),'fro')/norm(X_observed(Known_train),'fro');

    rmse_test_CGLS_rating = norm(X_CGLS_predict_theta(Known_test) - X_test(Known_test), 'fro') / sqrt(length(Known_test));
    rmse_train_CGLS_rating = norm(X_CGLS_predict_theta(Known_train) - X_observed(Known_train), 'fro') / sqrt(length(Known_train));

    rmse_test_CGLS_org = norm(Xcg2_theta(Known_test) - X_test(Known_test), 'fro') / sqrt(length(Known_test));
    rmse_train_CGLS_org = norm(Xcg2_theta(Known_train) - X_observed(Known_train), 'fro') / sqrt(length(Known_train));

    save_file = [result_root, '/CG_LS'];
    save(save_file, 'hist_CGLS_theta', 'out_time_CGLS', 'rte_test_CGLS', 'rte_train_CGLS', 'rmse_test_CGLS_rating', 'rmse_train_CGLS_rating', 'rmse_test_CGLS_org', 'rmse_train_CGLS_org')
end

%% THE END

if isempty(strcmp(computer, 'WIN'))
matlabpool('close')
end