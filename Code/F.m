function [f,f_principal] = F(prob,x,theta)
%F Cost function 
%

switch prob.loss_flag
        
    case 'LS'
        err = x.on_omega - prob.data; f_principal = 0.5*(err'*err);
        f = f_principal;
        
    case 'groupthre'
        delta = prob.delta;
        bridgingthre_ids = prob.bridgingthre_ids;
        
        f = loss_groupthre_matlab(x.on_omega,prob.data,theta,delta,bridgingthre_ids);
        
        f2 = loss_matlab(x.on_omega,prob.data,theta,delta);
        a = 1;
        
%         num_group = length(bridgingthre) + 1;
%         start_theta_id = 1;
%         for i_group=1:num_group
%             bridgingthre_id = bridgingthre_ids(i_group);
%             if i_group == 1
%                 sel_point_id = find(prob.data <= bridgingthre_id);
%             elseif i_group == num_group
%                 sel_point_id = find(prob.data > bridgingthre_id);
%             else
%                 sel_point_id = intersect(find(prob.theta <= bridgingthre_id), find(prob.theta > bridgingthre_ids(i_group-1)));
%             end
%             
%             f = f + loss_matlab(x.on_omega(sel_point_id),prob.data(sel_point_id),theta(start_theta_id:bridgingthre_id),delta);
%        
%         end
        
    case 'immthre'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        
        f = loss_immthre_matlab(x.on_omega,prob.data,theta,delta);
        
    case 'RMM_new'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        
        f = loss_matlab(x.on_omega,prob.data,theta,delta);
        
%         f_principal = 0;
%         for r=1:size(theta)
%            T = ones(size(x.on_omega));
%            T(prob.data>r) = -1;
%            err = max(0, delta - T .* (theta(r) - x.on_omega));
%            f_principal = f_principal + 0.5 * (err' * err);
%         end
%         
%         norm(f - f_principal,'fro')
%         a = 1;
        
    case 'RMM'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        f_principal = 0;
        for r=1:size(theta)
           T = ones(size(x.on_omega));
           T(prob.data>r) = -1;
%            err = max(0, 1 - T .* (theta(r) - x.on_omega));
           err = max(0, delta - T .* (theta(r) - x.on_omega));
           f_principal = f_principal + 0.5 * (err' * err);
        end
        f = f_principal;
    case 'RMM-revised'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        
        f_principal = 0;
        for r=1:size(theta)
           T = ones(size(x.on_omega));
           T(prob.data>prob.rating_levels(r)) = -1;
           err = max(0, delta - T .* (theta(r) - x.on_omega));
           f_principal = f_principal + 0.5 * (err' * err);
        end
        f = f_principal;
    case '1MM'
        err = max(0, 1 - x.on_omega .* prob.data);
        f = 0.5 * (err' * err);
end


if prob.mu>0
  f = f + 0.5*prob.mu*norm(1./x.sigma)^2 + 0.5*prob.mu*norm(x.sigma)^2;
end