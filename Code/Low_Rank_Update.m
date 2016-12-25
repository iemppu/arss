function [U S V] = Low_Rank_Update(U,S,V,U_L,S_L,V_L,k)
%%
% Low_Rank_Update for X = USV' + AB' = USV' + U_L*S_L*V_L'
%%
% R1 = R + a*b';
% [U1 S1 V1] = svds(R1,5);
s_l = diag(S_L);
for t = 1:k
    a = U_L(:,t);
    b = V_L(:,t)*s_l(t);
    m = U'*a;
    p = a-U*m;
    Ra = norm(p);
    P = 1/Ra*p;
    n = V'*b;
    q = b-V*n;
    Rb = norm(q);
    Q = 1/Rb*q;
    s_u = diag(S);
    s_u = [s_u;0];
    K = diag(s_u);
    K1 = K + [m;Ra]*[n;Rb]';
    [Uk,S,Vk] = svds(K1,length(s_u));
    U = [U P]*Uk;
    V = [V Q]*Vk;
end