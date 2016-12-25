function prob = make_prob_real_data(Data,Omega,n1,n2)
% L*R' is the exact data
% Omega is the sampling set; a sparse matrix.
% k is the rank to reconstruct


prob.n1 = n1; 
prob.n2 = n2; 

% prob.r = k;

prob.Omega = Omega;
[prob.Omega_i, prob.Omega_j] = ind2sub([prob.n1,prob.n2], Omega);
prob.m = length(Omega);
% noise = randn(prob.m,1);

% prob.data = Data(Omega); 
prob.data = Data;
% % noise_level = sigma*norm(prob.data)/norm(noise);
% % prob.data = prob.data + +noise_level*noise;
prob.temp_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data,prob.n1,prob.n2,prob.m);
  
% regularization parameter
%prob.mu = 1e-14;
prob.mu = 0;
