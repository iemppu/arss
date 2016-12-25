function g = grad(prob,x,theta)
%GRAD   Computes the gradient on the manifold
%
%  Computes the gradient at a point x of the cost function on 
%  the manifold. + 0.01*x.U +0.01*x.V

% multiple loss functions
switch prob.loss_flag
        
    case 'LS'
        d = x.on_omega - prob.data;
        
    case 'groupthre'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        bridgingthre_ids = prob.bridgingthre_ids;
        
        d = grad_groupthre_matlab(x.on_omega,prob.data,theta,delta,bridgingthre_ids);
        
        d2 = grad_matlab(x.on_omega,prob.data,theta,delta);
        a = 1;
        
    case 'immthre'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        
        d = grad_immthre_matlab(x.on_omega,prob.data,theta,delta);
        
    case 'RMM_new'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        
        d = grad_matlab(x.on_omega,prob.data,theta,delta);
        
%         d2 = zeros(size(x.on_omega));
%         for r=1:length(theta)
%            T = ones(size(x.on_omega));
%            T(prob.data>r) = -1;
%            err = T .* max(0, delta - T .* (theta(r) - x.on_omega));
%            d2 = d2 + err;
%         end
%         
%         norm(d2-d,'fro')
%         a = 1;
        
    case 'RMM'
        if nargin < 3
            theta = prob.theta;
        end
        delta = prob.delta;
        d = zeros(size(x.on_omega));
        for r=1:length(theta)
           T = ones(size(x.on_omega));
           T(prob.data>r) = -1;
%            err = T .* max(0, 1 - T .* (theta(r) - x.on_omega));
           err = T .* max(0, delta - T .* (theta(r) - x.on_omega));
           d = d + err;
        end
    case 'RMM-revised'
        if nargin < 3
            theta = prob.theta;
        end
        d = zeros(size(x.on_omega));
        for r=1:length(theta)
           T = ones(size(x.on_omega));
           T(prob.data>prob.rating_levels(r)) = -1;
           err = T .* max(0, 1 - T .* (theta(r) - x.on_omega));
           d = d + err;
        end
    case '1MM'
        err = max(0, 1 - x.on_omega .* prob.data);
        d = -prob.data .* err;
end

% d is the gradient in the Euclidean space
% below is the code to compute the gradient on the manifold

updateSval(prob.temp_omega, d, length(d));
prob.temp_omega = prob.temp_omega;

T = prob.temp_omega*x.V;
g.M = x.U'*T; 
g.Up = T - x.U*(x.U'*T);

g.Vp = prob.temp_omega'*x.U; 
g.Vp = g.Vp - x.V*(x.V'*g.Vp);

g.M = g.M + prob.mu*diag(x.sigma-1./x.sigma.^3);
    


