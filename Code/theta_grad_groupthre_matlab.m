function [theta_grad] = theta_grad_groupthre_matlab(X,Y,theta,delta,bridgingthre_ids)
% Compute thre gradient of f r.w.t. \theta

num_group = length(bridgingthre_ids) + 1;

index_bridgingthre = 1;


L = length(theta);
theta_grad = zeros(size(theta));
for i_theta = 1:L
    
    % check if current theta is a bridging threshold
    if i_theta == bridgingthre_ids(index_bridgingthre)
        if index_bridgingthre == 1 % if this is the first bridging threshold
            start_rating_value = 1;
            end_rating_value = bridgingthre_ids(index_bridgingthre);
        elseif index_bridgingthre == length(bridgingthre_ids) % if this is the last bridging threshold
            start_rating_value = bridgingthre_ids(index_bridgingthre-1);
            end_rating_value = length(theta) + 1;
        else % if this is a normal bridging threshold
            start_rating_value = bridgingthre_ids(index_bridgingthre-1);
            end_rating_value = bridgingthre_ids(index_bridgingthre+1);
        end
        
    
    else % if this is not a bridging threshold
        if i_theta < bridgingthre_ids(1) % if current threshold is in the first group
            start_rating_value = 1;
            end_rating_value = bridgingthre_ids(1);
        elseif i_theta > bridgingthre_ids(end); % if current threshold is in the last group
            start_rating_value = bridgingthre_ids(end) + 1;
            end_rating_value = length(theta) + 1;
        else % if current threshold is in the middle groups
            start_rating_value = bridgingthre_ids(index_bridgingthre-1);
            end_rating_value = bridgingthre_ids(index_bridgingthre);
        end
    end
    
    sel_point_id = intersect(find(Y>=start_rating_value),find(Y<=end_rating_value));
    
    if i_theta == bridgingthre_ids(index_bridgingthre)
        theta_grad(i_theta) = 1/2 * theta_grad_groupthre_sub_matlab(X(sel_point_id),Y(sel_point_id),theta,delta,i_theta);
        
        if index_bridgingthre < length(bridgingthre_ids)
            index_bridgingthre = index_bridgingthre + 1;
        end
    
    else
        theta_grad(i_theta) = theta_grad_groupthre_sub_matlab(X(sel_point_id),Y(sel_point_id),theta,delta,i_theta);
    end
    
end

