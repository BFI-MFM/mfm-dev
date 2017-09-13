function [elas1,elas2] = SEimfd1e(model,domain,bc,param,eig,X0,Nu)
% Calculating shock elasticity using implicit finite difference
% by Yiran Fan
%
% The input of the function:
%
% * model:  a structure that specifies the underlying process 
%                   d X_t = muX dt + sigmaX dW_t
%           the (log) cash flow process 
%                   d log C_t = muC dt + sigmaC dW_t
%           and the (log) SDF process 
%                   d log S_t = muS dt + sigmaS dW_t
%           Here muX, sigmaX, muC, sigmaC, muS, sigmaS are functions.
%           sigmaX, sigmaC, sigmaS are allowed to return row vectors.
%
%           But notice, if a vector is passed to sigmaX (etc), sigmaX
%           must output a matrix, where the first row is for the
%           corresponding results for the first element in the passed in
%           vector, so on so forth.
%
% * domain: a structure that specifies the domain of X and t, 
%           where X is a two-element vector gives the lower and upper 
%           bounds of X; while t is a multi-element vector gives all values
%           of t of interests; and nX gives the # of points of X.
%
% * bc:     a structure that specifies the boundary conditions for the lower
%           and upper bounds (bc.l and bc.u respectively). All boundary
%           conditions must in the format 
%                   0 = a0 + a1 phi(x) + a2 phi'(x) + a3 phi''(x),
%           where a0, a1, a2 can be functions of t. Therefore, for each
%           element in bc, say bc.l.C, there must be bc.l.C.const for
%           a0(t), bc.l.C.level for a1(t), bc.l.C.deriv1 for a2(t), 
%           and bc.l.C.deriv2 for a3(t).
%
% * param:  parameters in model
% 
% * eig:    eigenfunctions and eigenvalues for C and SC 
%           eig.C.e, eig.C.dloge, eig.C.eta
%           eig.SC.e, eig.SD.dloge, eig.SC.eta
%
% * X0:     a row vector
%
% * Nu:     a function that output a row vector for selecting shocks;
%           If input is a vector, then the function gives a matrix, where
%           the first row corresponds to the first element of the input
%           vector.

%% Step 1 Read the input

muX = model.muX; 
sigmaX = model.sigmaX;
betaC = model.muC;
alphaC = @(x) model.sigmaC(x).';
betaS = model.muS;
alphaS = @(x) model.sigmaS(x).';
betaSC = @(x) betaC(x) + betaS(x);
alphaSC = @(x) alphaC(x) + alphaS(x);
xx = linspace(domain.x(1), domain.x(end),domain.nx).';

dt = domain.dt;
tt = dt:dt:domain.T;
nx = domain.nx;
dx = xx(2) - xx(1);
sigmax = sigmaX(X0); 
nux = Nu(X0);
nX0 = size(X0,1); 

%% Step 2 Solve Shock Exposure Elasticity

alphax = alphaC(X0); % #shocks by #X0

aa = @(x) 1/2.*diag(sigmaX(x)*sigmaX(x).');
bb = @(x) muX(x)+diag(sigmaX(x)*alphaC(x)) + diag(sigmaX(x)*sigmaX(x).').*eig.C.dloge(x) ;

A = aa(xx);
B = bb(xx);
C = zeros(nx,1);

%% Step 2.1 First Type

disp('Calculating shock exposure elasticity for the first type');
tic;
psi0 = 1./eig.C.e(xx);

psi = CEimfd1(xx,tt,A,B,C,psi0,[X0;X0+dx;X0-dx],bc.C);
exp1 = psi  .*  repmat(exp(eig.C.eta*[0 tt]),3*nX0,1)    ...
        .*  repmat(eig.C.e([X0;X0+dx;X0-dx]),1,max(size(tt))+1);

% Shock elasiticy
tmp1 = exp1((nX0+1):(2*nX0),:);
tmp2 = exp1((2*nX0+1):(3*nX0),:);
exp1 = exp1(1:nX0,:);
tmpd = (tmp1-tmp2)/2/dx;
tmpld = tmpd ./ exp1; % #X0 by tt+1

nusigmax = diag( nux * sigmax.' ); %( #X0 by #shocks ) * (#X0 by #shock )'
nualphax = diag( nux * alphax ); % #X0 by 1
nusimgax = repmat(nusigmax,1,max(size(tt))+1); % #X0 by tt+1
nualphax = repmat(nualphax,1,max(size(tt))+1); % #X0 by tt+1

elas1.g = nusimgax.*tmpld+nualphax;
toc;   

%% Step 2.2 Second Type

disp('Calculating shock exposure elasticity for the second type');
tic;

psi0 = diag(Nu(xx) * alphaC(xx))./eig.C.e(xx);
psi = CEimfd1(xx,tt,A,B,C,psi0,X0,bc.C);
exp2 = psi  .*  repmat(exp(eig.C.eta*[0 tt]),nX0,1)    ...
        .*  repmat(eig.C.e(X0),1,max(size(tt))+1);

elas2.g = exp2 ./ exp1;
toc;

%% Step 3 Solve Shock Price Elasticity

alphax = alphaSC(X0); % column vector

aa = @(x) 1/2.*diag(sigmaX(x)*sigmaX(x).');
bb = @(x) muX(x)+diag(sigmaX(x)*alphaSC(x)) + diag(sigmaX(x)*sigmaX(x).').*eig.SC.dloge(x) ;

A = aa(xx);
B = bb(xx);
C = zeros(nx,1);

%% Step 3.1 First Type

disp('Calculating shock price elasticity for the first type');
tic;
psi0 = 1./eig.SC.e(xx);

psi = CEimfd1(xx,tt,A,B,C,psi0,[X0;X0+dx;X0-dx],bc.SC);
exp1 = psi  .*  repmat(exp(eig.SC.eta*[0 tt]),3*nX0,1)    ...
        .*  repmat(eig.SC.e([X0;X0+dx;X0-dx]),1,max(size(tt))+1);

% Shock elasiticy
tmp1 = exp1((nX0+1):(2*nX0),:);
tmp2 = exp1((2*nX0+1):(3*nX0),:);
exp1 = exp1(1:nX0,:);
tmpd = (tmp1-tmp2)/2/dx;
tmpld = tmpd ./ exp1; % #X0 by tt+1

nusigmax = diag( nux * sigmax.' ); %( #X0 by #shocks ) * (#X0 by #shock )'
nualphax = diag( nux * alphax ); % #X0 by 1
nusimgax = repmat(nusigmax,1,max(size(tt))+1); % #X0 by tt+1
nualphax = repmat(nualphax,1,max(size(tt))+1); % #X0 by tt+1

elas1.sg = nusimgax.*tmpld+nualphax;
toc;   

%% Step 3.2 Second Type

disp('Calculating shock price elasticity for the second type');
tic;

psi0 = diag(Nu(xx) * alphaSC(xx))./eig.SC.e(xx);
psi = CEimfd1(xx,tt,A,B,C,psi0,X0,bc.SC);
exp2 = psi  .*  repmat(exp(eig.SC.eta*[0 tt]),nX0,1)    ...
        .*  repmat(eig.SC.e(X0),1,max(size(tt))+1);

elas2.sg = exp2 ./ exp1;
toc;

%% Step 4 Calculate final elasticity

elas1.out = elas1.g - elas1.sg;
figure; plot([0 tt],elas1.g);title('First type Shock Exposure Elasticity')
figure; plot([0 tt],elas1.out);title('First type Shock Price Elasticity')
elas2.out = elas2.g - elas2.sg;
figure; plot([0 tt],elas2.g);title('Second type Shock Exposure Elasticity')
figure; plot([0 tt],elas2.out);title('Second type Shock Price Elasticity')