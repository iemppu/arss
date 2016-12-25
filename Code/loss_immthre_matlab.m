function [loss_val] = loss_immthre_matlab(X,Y,theta,delta)
% Compute the loss function: l = \sum_{ij} (\max(0,1-(\theta_{Y_{ij}-1} - X_{ij})) + \max(0,1-(\theta_{Y_{ij}} - X_{ij})))

L = length(theta);
f = 0;


% theta = [-inf;theta;inf];
theta = [theta(1)-100;theta;theta(end)+100];

l1 = delta - (theta(Y) - X);
l2 = delta - (X - theta(Y+1));

f = max(0,l1) + max(0,l2);


loss_val = f' * f * 0.5;

