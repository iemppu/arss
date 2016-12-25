function [theta] = symm_threshold(rating_levels,interval)
% Generate symmetric thresholds


num_theta = length(rating_levels) - 1;
start_value = -(num_theta / 2 * interval - 1/2 * interval);
end_value = -start_value;
theta = start_value:interval:end_value;



