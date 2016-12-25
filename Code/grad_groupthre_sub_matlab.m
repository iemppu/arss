function [grad] = grad_groupthre_sub_matlab(X,Y,theta,delta,theta_id)



grad = zeros(size(X));

for i_theta = theta_id
    T = ones(size(X));
    T(Y > i_theta) = -1;
    
    xi = max(0, delta - T .* (theta(i_theta) - X));
    
    grad = grad + xi .* T;
    
end


