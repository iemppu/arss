function prob = make_rating_prob_bias(L,R,Omega,k, sigma, rating_levels, bias)  % no scale and bias to make data
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
% prob.data = partXY(L', R',prob.Omega_i, prob.Omega_j,prob.m); 

% num_rand_perm = randperm(prob.m);
% prob.Omega_i = prob.Omega_i(num_rand_perm);
% prob.Omega_j = prob.Omega_j(num_rand_perm);
% prob.data  = prob.data(num_rand_perm);

noise_level = sigma*norm(prob.data)/norm(noise);
% prob.data = prob.data + noise_level*noise + bias;

num_levels = size(rating_levels, 1);

data_std = std(prob.data);
prob.data_std = data_std;
rho = linspace((2-num_levels)*data_std, (num_levels-2)*data_std, num_levels-1);

CC = 0.35;
rho = rho*CC;

max_val = max(prob.data);
min_val = min(prob.data);
interval_len = ( max_val - min_val ) / num_levels;

cuts_ = (min_val+interval_len):interval_len:(max_val-interval_len);
cuts = rho;
prob.Rdata = ones(size(prob.data));
for i=1:num_levels-1
    prob.Rdata(prob.data>cuts(i)) = i+1;
end

prob.cuts = cuts;

prob.raw_data = prob.data;
prob.data = prob.Rdata;

display('stop')
prob.temp_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data*1,prob.n1,prob.n2,prob.m);
  
% regularization parameter
%prob.mu = 1e-14;
prob.mu = 0;
