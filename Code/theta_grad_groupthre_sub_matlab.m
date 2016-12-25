function [theta_grad] = theta_grad_groupthre_sub_matlab(X,Y,theta,delta,i_theta)

L = length(theta);



% theta_grad = zeros(size(theta));






% for i_theta = rating_range
    T = ones(size(X));
    T(Y > i_theta) = -1;
    
    xi = max(0, delta - T .* (theta(i_theta) - X));
    
%     theta_grad(i_theta) = sum(xi.* (-T));
    theta_grad = sum(xi.* (-T));
% end

