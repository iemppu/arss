%%**************************************************
%% compute AT*z for 
%% A = beta1*X - beta2*Xold - beta3*Grad;
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun  
%%**************************************************

  function ATz = matvectransp_PP_my(z,param); 

  global Xh  Grad 

%   beta1 = param.beta1; 
%   beta2 = param.beta2; 
%   beta3 = param.beta3; 

  z1  = Xh.V*(Xh.U'*z); 
%   z2  = Xhold.V*(Xhold.U'*z); 
  ATz = (Grad'*z) - z1; 
%%**************************************************
