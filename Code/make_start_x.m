function x = make_start_x(prob)
% Compute a truncated SVD for starting guess.

lan_options.tol  = 1e-3;
M_omega = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data,prob.n1,prob.n2,prob.m);
[U,S,V] = lansvd(M_omega, prob.r,lan_options);
x.V = V;
x.sigma = diag(S);
x.U = U;

x = prepx(prob, x);


