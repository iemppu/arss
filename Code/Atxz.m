%%**************************************************
% compute A*z for 
% A = A_U * A_V' - Sparse_Z;
% 
% Atz = matvec(z,A_U,A_V,Sparse_Z); 
%
%%**************************************************
%%
  function Atz = Atxz(z); 

  global Sparse_Z
  
  Atz = (z' * Sparse_Z)';