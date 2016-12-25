function [x, opts]= make_start_x_adaptive(prob, opts)
% Compute a truncated SVD for starting guess.
global Grad Xh
lan_options.tol  = 1e-3;
Grad = sparse(prob.Omega_i, prob.Omega_j, ...
    prob.data,prob.n1,prob.n2,prob.m);

%%
% do patial svd using initialing guess prob.r = 2
prob.r = 2;
[U,S,V] = lansvd(Grad, prob.r,lan_options);

Xh.V = V;
Xh.sigma = diag(S);
Xh.U = U*diag(Xh.sigma);

sigular_value = diag(S);
sigular_max = sigular_value(1);
U_temp = U;
sigular_value_temp = sigular_value;

%
while (sigular_value(end)>opts.eta*sigular_max)
     prob.r = prob.r + 2;
     [U,S,V] = lansvd('matvec_PP_my','matvectransp_PP_my',size(U,1),size(V,1),2,'L',lan_options);
     Xh.V = [Xh.V V];
     Xh.sigma = diag(S);
     Xh.U = [Xh.U U*S];
     U_temp = [U_temp U];
     sigular_value = diag(S);
     sigular_value_temp = [sigular_value_temp; sigular_value];
     if prob.r > 0.5*min(size(Grad,1), size(Grad,2))
         break;
     end
end
Xh.V = Xh.V;
Xh.sigma = sigular_value_temp;
Xh.U = U_temp;
prob.r = length(sigular_value_temp);
opts.rank_increase = length(sigular_value_temp);
x = prepx(prob, Xh);


