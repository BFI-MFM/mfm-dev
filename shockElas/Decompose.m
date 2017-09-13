function [eta,e,xx] = Decompose(process,domain,bc,param)%,model)


%% Step 1 Read the input

muX = process.muX; 
sigmaX = process.sigmaX;
betaM = process.muM;
alphaM = @(x) process.sigmaM(x).';
xx = linspace(domain.x(1), domain.x(end),domain.nx).';

nx = domain.nx;
dx = xx(2) - xx(1);

%% Step 2 Create Generator

aa = @(x) 1/2.*diag(sigmaX(x)*sigmaX(x).');
bb = @(x) muX(x)+diag(sigmaX(x)*alphaM(x));
cc = @(x) betaM(x)+diag(alphaM(x).'*alphaM(x))./2;

A = aa(xx);
B = bb(xx);
C = cc(xx);

L = sparse(nx,nx);
for j = 2:(nx-1)
    
    L(j,j+1) = L(j,j+1) + A(j)/dx^2;
    L(j,j-1) = L(j,j-1) + A(j)/dx^2;
    L(j,j) = L(j,j) - 2*A(j)/dx^2;
    
    L(j,j+1) = L(j,j+1) + B(j)*(B(j)>0)/dx;
    L(j,j) = L(j,j) - B(j)*(B(j)>0)/dx;
    L(j,j-1) = L(j,j-1) - B(j)*(B(j)<0)/dx;
    L(j,j) = L(j,j) + B(j)*(B(j)<0)/dx;

%     L(j,j+1) = L(j,j+1) + B(j)/2/dx;
%     L(j,j-1) = L(j,j-1) - B(j)/2/dx;
%     
    L(j,j) = L(j,j) + C(j);
end
    
% Left boundary condition
L(1,:) = 0;
L(1,1) = L(1,1) + bc.l.level(0);
L(1,1) = L(1,1) - bc.l.deriv1(0)./dx;
L(1,2) = L(1,2) + bc.l.deriv1(0)./dx;
L(1,3) = L(1,3) + bc.l.deriv2(0)./dx^2;
L(1,1) = L(1,1) + bc.l.deriv2(0)./dx^2;
L(1,2) = L(1,2) - 2*bc.l.deriv2(0)./dx^2;
    
% Right boundary condition
L(nx,:) = 0;
L(nx,nx) = L(nx,nx) + bc.u.level(0);
L(nx,nx) = L(nx,nx) + bc.u.deriv1(0)./dx;
L(nx,nx-1) = L(nx,nx-1) - bc.u.deriv1(0)./dx;
L(nx,nx) = L(nx,nx) + bc.u.deriv2(0)./dx^2;
L(nx,nx-2) = L(nx,nx-2) + bc.u.deriv2(0)./dx^2;
L(nx,nx-1) = L(nx,nx-1) - 2*bc.u.deriv2(0)./dx^2;
    
%% Step 3 Decomposation

[cande,candeta] = eigs(L,10,0);

for j = 1:10
    numneg = sum(cande(:,j)<0);
    numpos = sum(cande(:,j)>0);
    if (numneg == nx)
        e = -cande(:,j);
        eta = candeta(j,j);
        disp(['Decomposation succeed, quit in ',num2str(j),' eigenvalue']);
        return;
    elseif (numpos == nx)
        e = cande(:,j);
        eta = candeta(j,j);
        disp(['Decomposation succeed, quit in ',num2str(j),' eigenvalue']);
        return;
    elseif (j == 10)
        dbstop if error;
        error('No all positive e');
    end
end   


