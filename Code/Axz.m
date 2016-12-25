%%**************************************************
% compute A*z for 
% A = A_U * A_V' - Sparse_Z;
% 
% Az = matvec(z,A_U,A_V,Sparse_Z); 
%
%%**************************************************
%%
  function Az = Axz(z); 

  global Sparse_Z
  
  Az =  Sparse_Z * z;