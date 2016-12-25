function [grad] = grad_matlab(X,Y,theta,delta)


L = length(theta);

grad = zeros(size(X));

for i_theta = 1:L
    T = ones(size(X));
    T(Y > i_theta) = -1;
    
    xi = max(0, delta - T .* (theta(i_theta) - X));
    
    grad = grad + xi .* T;
    
end


% for i_theta = 1:L
%     T = ones(size(X));
%     T(Y > i_theta) = -1;
%     
%     % indicator S
%     S = zeros(size(X));
%     S(delta >= T .* (theta(i_theta) - X)) = 1;
%     
%     xi = max(0, delta - T .* (theta(i_theta) - X));
%     
%     grad = grad + xi .* (T .* S);
%     
% end


