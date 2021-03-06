function [theta_grad] = theta_grad_matlab(X,Y,theta,delta)

L = length(theta);
theta_grad = zeros(size(theta));
for i_theta = 1:L
    T = ones(size(X));
    T(Y > i_theta) = -1;
    
    xi = max(0, delta - T .* (theta(i_theta) - X));
    
    theta_grad(i_theta) = sum(xi.* (-T));
end

