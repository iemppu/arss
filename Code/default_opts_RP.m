function opts = default_opts_RP()
% default options for LRGeomCG_pursuit
%
opts.maxit = 20;

% the ranks increased by each iteration
opts.rank_increase = 2;

% opts.rank_incr_search: options for the search of opts.rank_increase;
% if setting opts.rank_incr_search = 1, it will launch the rank search 
% process suggesedin the paper. This option is highly recommended to achieve 
% good performance
opts.rank_incr_search = 1;
opts.eta = 0.7;
% if setting opts.rank_incr_addap = 1, we adopt an addaptive updating of opts.rank_increase. 
opts.rank_incr_addap = 1;

%--------------------------------------------------------------------------
% stopping tollerance for the outer problem optimization
%--------------------------------------------------------------------------

% tollerance w.r.t. objective values
opts.rel_out_obj_tol = 0.01; 
% tollerance w.r.t. maximum singular value of the residual matrix
opts.rel_out_max_norm_tol = 0.001; 
% tollerance w.r.t. objective value difference
opts.rel_out_obj_diff_tol = 0.001; 

%--------------------------------------------------------------------------
% target options for MR 
%--------------------------------------------------------------------------

% Tolerance on the relative l_2 error on the sampling set Omega
opts.tar_maxit = 300;
opts.tar_rel_f_tol = 1e-6; % opts.tar_rel_f_tol = opts.rel_f_tol by default
opts.tar_rel_inner_tol = 0.001;  %  tollerance w.r.t. objective value difference

% Tolerance for detection of stagnation. 
opts.rel_tol_change_x = 1e-12;
opts.rel_tol_change_res = 1e-4;

% Verbosity 2 is very chatty.
opts.verbosity = 1;

% Strong Wolfe is needed in theory but not in practice and is a little
% slower
opts.strong_wolfe = false; 


%--------------------------------------------------------------------------
% stopping tollerances of LRGeomCG for master problem optimization (by Bart Vandereycken)
%--------------------------------------------------------------------------
opts.maxit = 30;
%  tollerance w.r.t. objective value difference
opts.rel_inner_tol = 0.01;  

% Tolerance on the Riemannian gradient of the objective function
opts.abs_grad_tol = 0; 
opts.rel_grad_tol = 1e-8;

% Tolerance on the l_2 error on the sampling set Omega
opts.abs_f_tol = 0;
opts.rel_f_tol = 1e-6;







