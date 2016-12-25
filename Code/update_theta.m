function [update_theta, update_f, theta_step_size] = update_theta(prob, x, theta, f_init, theta_size_init)
% update_theta
% For ratings, after updating matrix x, theta is required to be updated too
% 
% First compute the gradient for theta
% Then use the line search (Armijo search method) to set the step length
% Finally, return theta

grad_theta = compute_gradient_theta(prob, x, theta);
% dir = scaleTxM(grad_theta, -1);
dir = -grad_theta;

% tinit = exact_search_onlyTxM(prob, x, dir);
% tinit = 1/prob.m;

% tinit = theta_size_init/0.6;
tinit = theta_size_init;

% f_init = F(prob, x, theta);
[update_theta, update_f, succ, ~, iarm, theta_step_size] = armijo_search_for_theta(prob, x, theta, f_init, grad_theta, dir, tinit);

