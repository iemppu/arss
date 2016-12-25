function [loss_val] = loss_groupthre_sub_matlab(X,Y,theta,delta,theta_id)
% Compute the loss function: l = \sum_{ij} \sum_{z} (\xi_{ij}^{z})^{2}

f = 0;

for i_theta = theta_id
    
    % ---------- check if this threshold is bridging threshold ----------
    if i_theta == theta_id(1) || i_theta == theta_id(end) % if this is the first/last threshold in current group
        if i_theta ~= 1 && i_theta ~= length(theta) % if this is not the first/last threshold in all thresholds (then should be the bridging threshold)
            T = ones(size(X));
            T(Y > i_theta) = -1;
            err = max(0, delta - T .* (theta(i_theta) - X));
            f = f + 1/2 * (err' * err);
        end
    else % if this is a bridging threshold
        T = ones(size(X));
        T(Y > i_theta) = -1;
        err = max(0, delta - T .* (theta(i_theta) - X));
        f = f + (err' * err);
    end
    
end


loss_val = f * 0.5;

