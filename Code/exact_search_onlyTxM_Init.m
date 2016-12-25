function tmin = exact_search_onlyTxM_Init(prob,U,S,V, x, theta)

% Exact line search in the direction of dir on the tangent space of x
% !! so NOT the retracted curve, only use as guess !!

% x, current point

% dir is search direction
%
% returns: tmin


switch prob.loss_flag
    
    case 'LS'
        e_omega = x.on_omega - prob.data;
        
    case 'groupthre'
        if nargin < 4
            theta = prob.theta;
        end
        delta = prob.delta;
        
        bridgingthre_ids = prob.bridgingthre_ids;
        
        e_omega = grad_groupthre_matlab(x.on_omega,prob.data,theta,delta,bridgingthre_ids);
        
    case 'immthre'
        if nargin < 4
            theta = prob.theta;
        end
        delta = prob.delta;
        
        e_omega = grad_immthre_matlab(x.on_omega,prob.data,theta,delta);
        
    case 'RMM_new'
        if nargin < 4
            theta = prob.theta;
        end
        delta = prob.delta;
        
        e_omega = grad_matlab(x.on_omega,prob.data,theta,delta);
        
%         e_omega2 = zeros(size(x.on_omega));
%         for r=1:size(theta)
%            T = ones(size(x.on_omega));
%            T(prob.data>r) = -1;
%            err = T .* max(0, delta - T .* (theta(r) - x.on_omega));
%            e_omega2 = e_omega2 + err;
%         end
%         
%         norm(e_omega-e_omega2,'fro')
%         a = 1;
        
    case 'RMM'
        if nargin < 4
            theta = prob.theta;
        end
        delta = prob.delta;
        e_omega = zeros(size(x.on_omega));
        for r=1:size(theta)
           T = ones(size(x.on_omega));
           T(prob.data>r) = -1;
           err = T .* max(0, delta - T .* (theta(r) - x.on_omega));
           e_omega = e_omega + err;
        end
    case 'RMM-revised'
        if nargin < 4
            theta = prob.theta;
        end
        e_omega = zeros(size(x.on_omega));
        for r=1:size(theta)
           T = ones(size(x.on_omega));
           T(prob.data>prob.rating_levels(r)) = -1;
           err = T .* max(0, 1 - T .* (theta(r) - x.on_omega));
           e_omega = e_omega + err;
        end
    case '1MM'
        err = max(0, 1 - x.on_omega .* prob.data);
        e_omega = -prob.data .* err;
end

% beware: transposed!
dir_omega = partXY_blas(U',(V*S)', prob.Omega_i, prob.Omega_j, prob.m); 


% norm is f(t) = 0.5*||e+t*d||_F^2
% minimize analytically
% polynomial df/dt = a0+t*a1
% a0 = dir_omega*e_omega;
% a1 = dir_omega*dir_omega';
% tmin = -a0/a1;


A1 = dir_omega';
b = e_omega;
tmin = -2*A1\b;

