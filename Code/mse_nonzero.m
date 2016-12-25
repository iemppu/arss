% Mean square error   
% by yanyan
function [c] = mse_nonzero(a,b)
% a is recovered matrix
% b is testing matrix

  %c = full(sum(sum(abs(a-b).*(b>0)))./sum(sum(b>0)));
%   c = full(sum(sum((a.*(b>0)-b).^2))./sum(sum(b>0)));
c = full(sum(sum((a.*(b~=0)-b).^2))./sum(sum(b~=0)));
