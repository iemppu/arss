function prob = generate_prob(m,n,k,OS)
L0 = randn(m, k);
[L0,~] = qr(L0,0);
R = randn(n, k); 
[R,~] = qr(R,0);
S = diag(1000*(abs(randn(k,1)))); 
L = L0*S;

dof = k*(m+n-k);
rank_reconstruct = k;

% make random sampling, problem and initial guess
samples = floor(OS * dof);
Omega = make_rand_Omega(m,n,samples);
prob = make_prob(L,R,Omega,rank_reconstruct,0); % <- you can choose another rank here