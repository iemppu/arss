function [loss_val] = loss_matlab(X,Y,theta,delta)
% Compute the loss function: l = \sum_{ij} \sum_{z} (\xi_{ij}^{z})^{2}

L = length(theta);
f = 0;
for i_theta = 1:L
    T = ones(size(X));
    T(Y > i_theta) = -1;
    err = max(0, delta - T .* (theta(i_theta) - X));
    f = f + (err' * err);
end


loss_val = f * 0.5;

