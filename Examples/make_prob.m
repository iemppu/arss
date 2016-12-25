function prob = make_prob(L,R,Omega,k, sigma)
% L*R' is the exact data
% Omega is the sampling set; a sparse matrix.
% k is the rank to reconstruct


prob.n1 = size(L,1); 
prob.n2 = size(R,1); 
prob.r = k;

prob.Omega = Omega;
[prob.Omega_i, prob.Omega_j] = ind2sub([prob.n1,prob.n2], Omega);
prob.m = length(Omega);
noise = randn(prob.m,1);


prob.data = partXY(L', R',prob.Omega_i, prob.Omega_j,prob.m)'; 

% num_rand_perm = randperm(prob.m);
% prob.Omega_i = prob.Omega_i(num_rand_perm);
% prob.Omega_j = prob.Omega_j(num_rand_perm);
% prob.data  = prob.data(num_rand_perm);

noise_level = sigma*norm(prob.data)/norm(noise);
prob.data = prob.data +noise_level*noise;
prob.temp_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data*1,prob.n1,prob.n2,prob.m);
  
% regularization parameter
%prob.mu = 1e-14;
prob.mu = 0;
