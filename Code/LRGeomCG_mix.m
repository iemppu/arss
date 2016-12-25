function [x, histout, theta, itc, fail] = LRGeomCG_mix(prob, opts, x0, fc0, rel_inner_tol)
% LRGEOMCG    Low-rank matrix completion by Geometric CG. 
%   Uses Armijo rule, polynomial linesearch on embedded submanifold of
%   fixed rank matrices.
%           
% Input: prob     = problem instance, see MAKE_PROB.
%        opts     = options, see DEFAULT_OPTS.           
%        x0       = starting guess.
%
% Output: x       = solution.
%         histout = iteration history. Each row of histout is
%                   [rel_norm(grad), rel_err_on_omega, relative_change, ...
%                        number of step length reductions, restarts]
%         fail    = flag for failure
%
% See also SMALL_EXAMPLE, MAKE_PROB, DEFAULT_OPTS. 

% (C) Bart Vandereycken, 2011-2012.
% Adapted from steep by C. T. Kelley, Dec 20, 1996.
%
% Paper: Low-rank matrix completion by Riemannian optimization.
%
% Info and bugs: bart.vandereycken@epfl.ch

if isfield(prob,'theta')
    theta = prob.theta;
end

if nargin<5
    rel_inner_tol = 1e-10;
end
if nargin<4
    fc0 = F(prob,x0);
end

%beta_type = 'F-R';
beta_type = 'P-R';

fail = true;

norm_M_Omega = prob.norm_M_Omega;

ORTH_VALUE = 0.1; % the search directions should be almost orthogonal

itc = 1;
xc = x0;
fc = fc0;
gc = grad(prob,xc);
ip_gc = ip(xc,gc,gc);
fold = 2*fc; reschg = -1; 
% first search-dir is steepest gradient
dir = scaleTxM(gc,-1);
rel_grad = sqrt(ip_gc)/max(1,norm(xc.sigma));

ithist=zeros(opts.maxit,5);
beta = 0;

% COMPUTE RMSE
if prob.rmse_flag == 1;
    rating_levels = prob.rating_levels;
    % train_rmse_wrong = norm(xc.on_omega-prob.data, 'fro')/sqrt(prob.m);
    x_train = sparse(prob.Omega_i, prob.Omega_j, prob.data, prob.n1, prob.n2);

    x_predict_rowMatrix = xc.U * diag(xc.sigma) * xc.V';

    % discrete prediction
    if strcmp(prob.loss_flag, 'LS')
        x_predict = x_predict_rowMatrix;
            x_predict(x_predict_rowMatrix>=0) = 1;
            x_predict(x_predict_rowMatrix<0) = -1;
    else
    cuts = theta;
    x_predict = ones(size(x_predict_rowMatrix));
    x_predict = x_predict * rating_levels(1);
    for i=1:size(rating_levels,1)-1
    %     x_predict(x_predict_rowMatrix>cuts(i)) = i+1;
            x_predict(x_predict_rowMatrix>cuts(i)) = rating_levels(i+1);
    end
    end
    train_rmse = sqrt(mse_nonzero(x_predict, x_train));
    ithist(itc,6) = train_rmse;
    train_rmse2 = norm(x_predict(prob.train_indices) - x_train(prob.train_indices), 'fro') / sqrt(prob.m);
    ithist(itc,7) = train_rmse2;

    test_rmse = sqrt(mse_nonzero(x_predict, prob.X_test_set));
    ithist(itc,8) = test_rmse;
    test_rmse2 = norm(x_predict(prob.test_indices) - prob.X_test_set(prob.test_indices), 'fro') / sqrt(length(prob.test_indices));
    ithist(itc,9) = test_rmse2;
end

ithist(itc,1) = rel_grad;
ithist(itc,2) = (2*fc);%/norm_M_Omega
ithist(itc,3) = reschg;
ithist(itc,4) = 0;    
ithist(itc,5) = etime(clock,prob.t0);    

t_theta_init = 1/prob.m;

for itc=2:opts.maxit
  tinit = exact_search_onlyTxM(prob, xc,dir,theta);
  
%   fprintf('Iteration %d\n',itc)
  fc_b = fc;
  [xc_new,fc,succ,~,iarm] = armijo_search(prob, xc, fc, gc, dir, tinit, theta);
  
  if ~succ && beta ~= 0
    if opts.verbosity > 0; disp('Line search failed on CG. Resetting to gradient.'); end
    beta = 0;
    % if we did cg, reset to steepest descent and try again
    dir = scaleTxM(gc,-1);
%     tinit(itc) = exact_search_onlyTxM(prob, xc,dir,theta);
    tinit = exact_search_onlyTxM(prob, xc,dir,theta);
    [xc_new,fc,succ,~,iarm] = armijo_search(prob, xc,fc_b,gc,dir, tinit, theta);
  end
  
  %   update theta
% ========================================================================
  if isfield(prob,'theta') && isfield(prob,'theta_flag')
      if strcmp(prob.theta_flag, 'update') && ~strcmp(prob.loss_flag,'LS')
          fc_theta_init = fc;
          [new_theta, ft_theta, t_theta] = update_theta(prob, xc_new, theta, fc_theta_init, t_theta_init);
          theta = new_theta;
          prob.theta = theta;  % NOTICE: update theta to prob
          t_theta_init = t_theta;
          fc = ft_theta;
      end
  end
  
  % if it still fails (beta is always 0 here -> steepest descent)
  % then we give up
  if ~succ
    x = xc_new;
    

    % COMPUTE RMSE
    if prob.rmse_flag == 1;
        rating_levels = prob.rating_levels;
        % train_rmse_wrong = norm(xc.on_omega-prob.data, 'fro')/sqrt(prob.m);
        x_train = sparse(prob.Omega_i, prob.Omega_j, prob.data, prob.n1, prob.n2);

        x_predict_rowMatrix = xc.U * diag(xc.sigma) * xc.V';

        % discrete prediction
        if strcmp(prob.loss_flag, 'LS')
            x_predict = x_predict_rowMatrix;
                x_predict(x_predict_rowMatrix>=0) = 1;
                x_predict(x_predict_rowMatrix<0) = -1;
        else
        cuts = theta;
        x_predict = ones(size(x_predict_rowMatrix));
        x_predict = x_predict * rating_levels(1);
        for i=1:size(rating_levels,1)-1
        %     x_predict(x_predict_rowMatrix>cuts(i)) = i+1;
                x_predict(x_predict_rowMatrix>cuts(i)) = rating_levels(i+1);
        end
        end
        train_rmse = sqrt(mse_nonzero(x_predict, x_train));
        ithist(itc,6) = train_rmse;
        train_rmse2 = norm(x_predict(prob.train_indices) - x_train(prob.train_indices), 'fro') / sqrt(prob.m);
        ithist(itc,7) = train_rmse2;

        test_rmse = sqrt(mse_nonzero(x_predict, prob.X_test_set));
        ithist(itc,8) = test_rmse;
        test_rmse2 = norm(x_predict(prob.test_indices) - prob.X_test_set(prob.test_indices), 'fro') / sqrt(length(prob.test_indices));
        ithist(itc,9) = test_rmse2;
    end
    
    ithist(itc,1) = rel_grad;
    ithist(itc,2) = (2*fc);%/norm_M_Omega
    ithist(itc,3) = reschg;
    ithist(itc,4) = iarm;  
    ithist(itc,5) = etime(clock,prob.t0);  
    histout=ithist(1:itc,:); 
    if opts.verbosity > 0; disp('Line search failed on steepest descent. Exiting...'); end
    return
  end
  
  
  % grad(new x)
  gc_new = grad(prob,xc_new);
  ip_gc_new = ip(xc_new,gc_new,gc_new);
  

  % Test for stopping
  if (itc>3)&&((ithist(itc-2,2)-ithist(itc-1,2))/ithist(itc-2,2) < rel_inner_tol)
     break;
  end
  if sqrt(2*fc) < opts.abs_f_tol
    if opts.verbosity > 0
      disp('Abs f tol reached.')
    end
    fail = false;
    break;
  end
  if sqrt(2*fc)/norm_M_Omega < opts.rel_f_tol
    if opts.verbosity > 0; disp('Relative f tol reached.'); end
    fail = false;
    break;
  end
  
  rel_grad = sqrt(ip_gc_new)/max(1,norm(xc_new.sigma));
  if rel_grad < opts.rel_grad_tol
    if opts.verbosity > 0; disp('Relative gradient tol reached.'); end
    fail = false;
    break;
  end
  
  % for 'stagnation stopping criterion' after 10 iters
  reschg = abs(1-sqrt(2*fc)/sqrt(2*fold) );  % LMARank's detection
  %reschg = abs(sqrt(2*fc) - sqrt(2*fold)) / max(1,norm_M_Omega);
  if itc > 10 && reschg < opts.rel_tol_change_res
    if opts.verbosity > 0; disp('Iteration stagnated rel_tol_change_res.'); end
    fail = true;
    break;
  end
    
  % new dir = -grad(new x)
  %           + beta * vectTransp(old x, old dir, tmin*old dir)
  gc_old = transpVect(prob,xc,gc,xc_new,1);  
  dir = transpVect(prob,xc,dir,xc_new,1);
      
  % we test how orthogonal the previous gradient is with
  % the current gradient, and possibly reset the to the gradient
  ip_gc_old_new = ip(xc_new,gc_old,gc_new);
  orth_grads = ip_gc_old_new / ip_gc_new;
  
  if abs(orth_grads) >= ORTH_VALUE
    if opts.verbosity > 1; disp('New gradient is almost orthogonal to current gradient. This is good, so we can reset to steepest descent.'); end
    beta = 0;
    dir = plusTxM(gc_new, dir, -1, beta);
    
  else % Compute the CG modification
    % beta
    if strcmp(beta_type, 'F-R')  % Fletcher-Reeves
      beta = ip_gc_new / ip_gc;
      % new dir
      dir = plusTxM(gc_new, dir, -1, beta);
    elseif strcmp(beta_type, 'P-R')  % Polak-Ribiere+
      % vector grad(new) - transported grad(current)
      beta = (ip_gc_new - ip_gc_old_new) / ip_gc;            
      beta = max(0,beta);
      dir = plusTxM(gc_new, dir, -1, beta);
    end
  end
  
  % check if dir is descent, if not take -gradient (i.e. beta=0)
  g0 = ip(xc_new,gc_new,dir);
  if g0>=0
    if opts.verbosity > 1; 
      disp('New search direction not a descent direction. Resetting to -grad.');       
    end
    dir = scaleTxM(gc_new, -1);
    beta = 0;
  end
  
  
  % update _new to current
  gc = gc_new;
  ip_gc = ip_gc_new;
  xc = xc_new;
  fold = fc;
  

    % COMPUTE RMSE
    if prob.rmse_flag == 1;
        rating_levels = prob.rating_levels;
        % train_rmse_wrong = norm(xc.on_omega-prob.data, 'fro')/sqrt(prob.m);
        x_train = sparse(prob.Omega_i, prob.Omega_j, prob.data, prob.n1, prob.n2);

        x_predict_rowMatrix = xc.U * diag(xc.sigma) * xc.V';

        % discrete prediction
        if strcmp(prob.loss_flag, 'LS')
            x_predict = x_predict_rowMatrix;
                x_predict(x_predict_rowMatrix>=0) = 1;
                x_predict(x_predict_rowMatrix<0) = -1;
        else
        cuts = theta;
        x_predict = ones(size(x_predict_rowMatrix));
        x_predict = x_predict * rating_levels(1);
        for i=1:size(rating_levels,1)-1
        %     x_predict(x_predict_rowMatrix>cuts(i)) = i+1;
                x_predict(x_predict_rowMatrix>cuts(i)) = rating_levels(i+1);
        end
        end
        train_rmse = sqrt(mse_nonzero(x_predict, x_train));
        ithist(itc,6) = train_rmse;
        train_rmse2 = norm(x_predict(prob.train_indices) - x_train(prob.train_indices), 'fro') / sqrt(prob.m);
        ithist(itc,7) = train_rmse2;

        test_rmse = sqrt(mse_nonzero(x_predict, prob.X_test_set));
        ithist(itc,8) = test_rmse;
        test_rmse2 = norm(x_predict(prob.test_indices) - prob.X_test_set(prob.test_indices), 'fro') / sqrt(length(prob.test_indices));
        ithist(itc,9) = test_rmse2;
    end
    
  ithist(itc,1) = rel_grad;
  ithist(itc,2) = (2*fc);
  ithist(itc,3) = reschg;
  ithist(itc,4) = iarm;    
  ithist(itc,5) = etime(clock,prob.t0);  
end

%/norm_M_Omega
x = xc_new;

% COMPUTE RMSE
if prob.rmse_flag == 1;
    rating_levels = prob.rating_levels;
    % train_rmse_wrong = norm(xc.on_omega-prob.data, 'fro')/sqrt(prob.m);
    x_train = sparse(prob.Omega_i, prob.Omega_j, prob.data, prob.n1, prob.n2);

    x_predict_rowMatrix = xc.U * diag(xc.sigma) * xc.V';

    % discrete prediction
    if strcmp(prob.loss_flag, 'LS')
        x_predict = x_predict_rowMatrix;
            x_predict(x_predict_rowMatrix>=0) = 1;
            x_predict(x_predict_rowMatrix<0) = -1;
    else
    cuts = theta;
    x_predict = ones(size(x_predict_rowMatrix));
    x_predict = x_predict * rating_levels(1);
    for i=1:size(rating_levels,1)-1
    %     x_predict(x_predict_rowMatrix>cuts(i)) = i+1;
            x_predict(x_predict_rowMatrix>cuts(i)) = rating_levels(i+1);
    end
    end
    train_rmse = sqrt(mse_nonzero(x_predict, x_train));
    ithist(itc,6) = train_rmse;
    train_rmse2 = norm(x_predict(prob.train_indices) - x_train(prob.train_indices), 'fro') / sqrt(prob.m);
    ithist(itc,7) = train_rmse2;

    test_rmse = sqrt(mse_nonzero(x_predict, prob.X_test_set));
    ithist(itc,8) = test_rmse;
    test_rmse2 = norm(x_predict(prob.test_indices) - prob.X_test_set(prob.test_indices), 'fro') / sqrt(length(prob.test_indices));
    ithist(itc,9) = test_rmse2;
end

ithist(itc,1) = rel_grad;
ithist(itc,2) = (2*fc);
ithist(itc,3) = reschg;
ithist(itc,4) = iarm; 
ithist(itc,5) = etime(clock,prob.t0);  
histout=ithist(1:itc,:); 

