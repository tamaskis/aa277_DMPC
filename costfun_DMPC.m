function J = costfun_DMPC(f,xk1,xk2,Q,R,P,m,u1,u2,dt)
%Inputs - 
%   f: dynamics function
%   xk1: state of spacecraft 1 at time k [a; ex; ey; i; omega; u;]
%   xk2: state of spacecraft 2 at time k
%   Q,R,P weights on state, input, and terminal state respectively
%       Q  = [q1,...,q6], R = [r1,...,r3], P = [p1,...,p6]
%
%   m: prediction horizon (number of steps)
%   u1: vector of control inputs for spacecraft 1
%   u2: vector of control inputs for spacecraft 2
%       *when implementing, one of u1,u2 will be the optimal control input
%       sequence from the previous time step
%   dt: time step

%check that u and x have the correct number of elements
if (size(u1,2) ~= m) || (size(u2,2) ~= m)
    error('length(u,x) ~= m');
end

%initialize predicted states
x1 = xk1(:,1); 
x2 = xk2(:,1);

%propogate through the horizon
for i = 1:m-1
    x1(:,i+1) = f(x1(:,i),u1(:,i),dt);
    x2(:,i+1) = f(x2(:,i),u2(:,i),dt);
end

%compute cost (summed up over prediction horizon)
J1 = sum(Q*(x1(:,1:end-1)-x2(:,1:end-1)).^2 + R*u1(:,1:end-1).^2) + P*(x1(:,end)-x2(:,end)).^2;
J2 = sum(Q*(x1(:,1:end-1)-x2(:,1:end-1)).^2 + R*u2(:,1:end-1).^2) + P*(x1(:,end)-x2(:,end)).^2;
J = J1+J2;
end
