function [grad] = grad_immthre_matlab(X,Y,theta,delta)
% Compute the gradient of immediate threshold loss w.r.t. X

L = length(theta);

% grad = zeros(size(X));

% theta = [-inf;theta;inf];
theta = [theta(1)-100;theta;theta(end)+100];

l1 = delta - (theta(Y) - X);
l2 = delta - (X - theta(Y+1));

xi = max(0,l1) + max(0,l2);

T1 = zeros(size(xi));
T1(X >= theta(Y)-delta) = 1;

T2 = zeros(size(xi));
T2(X <= theta(Y+1)+delta) = -1;

T = T1 + T2;

grad = T .* xi;


