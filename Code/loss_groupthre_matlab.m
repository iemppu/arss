function [loss_val] = loss_groupthre_matlab(X,Y,theta,delta,bridgingthre_ids)
% Compute the loss function: l = \sum_{ij} (\max(0,1-(\theta_{Y_{ij}-1} - X_{ij})) + \max(0,1-(\theta_{Y_{ij}} - X_{ij})))

loss_val = 0;

num_group = length(bridgingthre_ids) + 1;
start_theta_id = 1;
for i_group=1:num_group
    if i_group == 1
        bridgingthre_id = bridgingthre_ids(i_group);
        sel_point_id = find(Y <= bridgingthre_id);
        end_theta_id = bridgingthre_id;
    elseif i_group == num_group
        sel_point_id = find(Y > bridgingthre_id);
        end_theta_id = length(theta);
    else
        bridgingthre_id = bridgingthre_ids(i_group);
        sel_point_id = intersect(find(Y <= bridgingthre_id), find(Y > bridgingthre_ids(i_group-1)));
        end_theta_id = bridgingthre_ids(i_group);
    end

%     loss_val = loss_val + loss_matlab(X(sel_point_id),Y(sel_point_id),theta(start_theta_id:bridgingthre_id),delta);
%     loss_val = loss_val + loss_matlab(X(sel_point_id),Y(sel_point_id),theta(start_theta_id:end_theta_id),delta);
    loss_val = loss_val + loss_groupthre_sub_matlab(X(sel_point_id),Y(sel_point_id),theta,delta,(start_theta_id:end_theta_id));
    start_theta_id = bridgingthre_id;

end

