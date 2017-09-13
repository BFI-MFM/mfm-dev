function [tt ts] = TSimfd1(model,domain,bc,param,X0,modelname)

%% Step 1 Read the input

muX = model.muX; 
sigmaX = model.sigmaX;
betaS = model.muS;
alphaS = @(x) model.sigmaS(x).';

xx = linspace(domain.x(1), domain.x(end),domain.nx).';
dt = domain.dt;
tt = dt:dt:domain.T;
nx = domain.nx;
dx = xx(2) - xx(1);

%% Step 2 Solve Expectation

aa = @(x) 1/2.*diag(sigmaX(x)*sigmaX(x).');
bb = @(x) muX(x)+diag(sigmaX(x)*alphaS(x));
cc = @(x) betaS(x)+diag(alphaS(x).'*alphaS(x))./2;

A = aa(xx);
B = bb(xx);
C = cc(xx);

phi0 = ones(nx,1);
expectation = CEimfd1(xx,tt,A,B,C,phi0,X0,bc);

ts = -1./repmat(tt,size(X0,1),1) .* log(expectation(:,2:end));
