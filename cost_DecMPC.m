function J = cost_DecMPC( xk, xref, Q, R, P, N, u, dt )
%Inputs - 
%   f: dynamics function
%   xk: state at time k [a; ex; ey; i; omega; u;
%   xref: reference trajectory
%   Q,R,P weights on state, input, and terminal state respectively
%       Q  = [q1,...,q6], R = [r1,...,r3], P = [p1,...,p6]
%
%   m: prediction horizon (number of steps)
%   u: vector of control inputs
%   dt: time step

%check that u and xref have the correct number of elements
if (length(u) ~= N) || (length(xref) ~= N)
    error('length(u) ~= m');
end

x = xk(:,1); 
%populate prediction horizon
for i = 1:N-1
    x(:,i+1) = nonlinear_dynamics(x(:,i),u(:,i),dt);
end

%compute cost (summed up over prediction horizon)
J = sum(Q*(x(:,1:end-1)-xref(:,1:end-1)).^2 + ...
        R*u(:,1:end-1).^2) + ...
		P*(x(:,end)-xref(:,end)).^2;
end