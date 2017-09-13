function expectation = CEimfd1(xx,tt,A,B,C,phi0,X0,bc)
% This function solves the following PDE using implicit finite difference:
%           d/dt phi = C phi + B d/dx phi + A d^2/dx^2 phi
%
% Inputs are:
% * xx:     The grid of X, must be COLUMN vector
% * tt:     The grid of T;
% * A,B,C:  Corresponding coefficients in the PDE, each should be a column 
%           vector of the same lenght of xx. These are not functions of t.
% * phi0:   The initial condition, should be a column 
%           vector of the same length of xx.
% * X0:     The points wanted to monitor, should be a COLUMN vector
% * bc:     The boundary condition of the following form:
%                   0 = a0 + a1 phi(x) + a2 phi'(x) + a3 phi''(x),
%           where a0, a1, a2 can be functions of t.

%% Section 1 Read Input

dx = xx(2) - xx(1);
dt = tt(2) - tt(1);
nx = max(size(xx));

%% Section 2 Construct operator

L = sparse(-eye(nx));
for j = 2:(nx-1)
    
    L(j,j+1) = L(j,j+1) + dt*A(j)/dx^2;
    L(j,j-1) = L(j,j-1) + dt*A(j)/dx^2;
    L(j,j) = L(j,j) - dt*2*A(j)/dx^2;
    
    L(j,j+1) = L(j,j+1) + dt*B(j)*(B(j)>0)/dx;
    L(j,j) = L(j,j) - dt*B(j)*(B(j)>0)/dx;
    L(j,j-1) = L(j,j-1) - dt*B(j)*(B(j)<0)/dx;
    L(j,j) = L(j,j) + dt*B(j)*(B(j)<0)/dx;

%     L(j,j+1) = L(j,j+1) + dt*B(j)/2/dx;
%     L(j,j-1) = L(j,j-1) - dt*B(j)/2/dx;
    
    L(j,j) = L(j,j) + dt*C(j);
end

%% Section 3 Solve PDE

% t = 0 case
phi = phi0;
expectation = spline(xx(~isnan(phi0)),phi0(~isnan(phi0)),X0);

for t = tt
    b = -phi;
    
    % Left boundary condition
    L(1,:) = 0;
    L(1,1) = L(1,1) + bc.l.level(t);
    L(1,1) = L(1,1) - bc.l.deriv1(t)./dx;
    L(1,2) = L(1,2) + bc.l.deriv1(t)./dx;
    L(1,3) = L(1,3) + bc.l.deriv2(t)./dx^2;
    L(1,1) = L(1,1) + bc.l.deriv2(t)./dx^2;
    L(1,2) = L(1,2) - 2*bc.l.deriv2(t)./dx^2;
    b(1) = - bc.l.const(t);
    
    % Right boundary condition
    L(nx,:) = 0;
    L(nx,nx) = L(nx,nx) + bc.u.level(t);
    L(nx,nx) = L(nx,nx) + bc.u.deriv1(t)./dx;
    L(nx,nx-1) = L(nx,nx-1) - bc.u.deriv1(t)./dx;
    L(nx,nx) = L(nx,nx) + bc.u.deriv2(t)./dx^2;
    L(nx,nx-2) = L(nx,nx-2) + bc.u.deriv2(t)./dx^2;
    L(nx,nx-1) = L(nx,nx-1) - 2*bc.u.deriv2(t)./dx^2;
    b(nx) = - bc.u.const(t);
    
    phi = L\b;
    tmp = spline(xx,phi,X0);
    expectation = [expectation tmp];
end

