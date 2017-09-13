function [tt,elas1,elas2] = SEimfd1(model,domain,bc,param,X0,Nu,modelname)
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
%           sigmaX, sigmaC, sigmaS are allowed to return row vectors, 
%           if a single value is passed in.
%
%           But notice, if a vector (must be column) is passed to sigmaX 
%           (etc), sigmaX must output a matrix, where the first row is for 
%           the corresponding results for the first element in the passed in
%           vector, so on so forth.
%
% * domain: a structure that specifies the domain of X and t. Having four
%           fields: 
%           x -- a two-element vector gives the lower and upper bounds of X; 
%           nx -- the number of grids to discretize X;
%           dt -- The length of time stamp;
%           T -- the termination time.
%
% * bc:     a structure that specifies the boundary conditions for the lower
%           and upper bounds (bc.l and bc.u respectively). All boundary
%           conditions must in the format 
%                   0 = a0 + a1 phi(x) + a2 phi'(x) + a3 phi''(x),
%           where a0, a1, a2 can be functions of t. Therefore, for each
%           element in bc, say bc.C.l, there must be bc.C.l.const for
%           a0(t), bc.C.l.level for a1(t), bc.C.l.deriv1 for a2(t), 
%           and bc.C.l.deriv2 for a3(t).
%
% * param:  parameters in model
% 
% * X0:     a column vector
%
% * Nu:     a function that output a row vector for selecting shocks;
%           If input is a vector, then the function gives a matrix, where
%           the first row corresponds to the first element of the input
%           vector.
% 
% * modelname: A string provides a directory where the automatically generated elasticity
%              plots should be saved to. If empty, the plots will not be saved.

%% Step 1 Read the input

muX = model.muX; % #X by 1
sigmaX = model.sigmaX; % #X by #shocks
betaC = model.muC; % #X by 1
alphaC = @(x) model.sigmaC(x).'; % #shocks by #X
betaS = model.muS; % #X by 1
alphaS = @(x) model.sigmaS(x).'; % #shocks by #X
betaSC = @(x) betaC(x) + betaS(x); % #X by 1
alphaSC = @(x) alphaC(x) + alphaS(x); % #shocks by #X
xx = linspace(domain.x(1), domain.x(end),domain.nx).'; 

if (~iscolumn(X0))
    X0 = X0.';
end    

dt = domain.dt;
tt = dt:dt:domain.T;
nx = domain.nx;
dx = xx(2) - xx(1);
sigmax = sigmaX(X0); % #X0 by #shocks
nux = Nu(X0); % #X0 by #shocks 
nX0 = size(X0,1); 

%% Step 2 Solve Shock Exposure Elasticity

alphax = alphaC(X0); % #shocks by #X0

aa = @(x) 1/2.*diag(sigmaX(x)*sigmaX(x).');
bb = @(x) muX(x)+diag(sigmaX(x)*alphaC(x));
cc = @(x) betaC(x)+diag(alphaC(x).'*alphaC(x))./2 ;

A = aa(xx);
B = bb(xx);
C = cc(xx);

%% Step 2.1 First Type

disp('Calculating shock exposure elasticity for the first type');
tic;
phi0 = ones(nx,1);

exp1 = CEimfd1(xx,tt,A,B,C,phi0,[X0;X0+dx;X0-dx],bc.C);

% Shock elasiticy
tmp1 = exp1((nX0+1):(2*nX0),:);
tmp2 = exp1((2*nX0+1):(3*nX0),:);
exp1 = exp1(1:nX0,:);
tmpd = (tmp1-tmp2)/2/dx;
tmpld = tmpd ./ exp1; % #X0 by tt+1

nusigmax = diag( nux * sigmax.' ); %( #X0 by #shocks ) * (#X0 by #shock )'
nualphax = diag( nux * alphax ); % #X0 by 1
nusigmax = repmat(nusigmax,1,max(size(tt))+1); % #X0 by tt+1
nualphax = repmat(nualphax,1,max(size(tt))+1); % #X0 by tt+1

elas1.g = nusigmax.*tmpld+nualphax;
toc;    

%% Step 2.2 Second Type

disp('Calculating shock exposure elasticity for the second type');
tic;

phi0 = diag(Nu(xx) * alphaC(xx));
exp2 = CEimfd1(xx,tt,A,B,C,phi0,X0,bc.C);
elas2.g = exp2 ./ exp1;

toc;

%% Step 3 Solve Shock Price Elasticity

alphax = alphaSC(X0);  % #shocks by #X0

aa = @(x) 1/2.*diag(sigmaX(x)*sigmaX(x).');
bb = @(x) muX(x)+diag(sigmaX(x)*alphaSC(x));
cc = @(x) betaSC(x)+diag(alphaSC(x).'*alphaSC(x))./2;

A = aa(xx);
B = bb(xx);
C = cc(xx);

%% Step 3.1 First Type

disp('Calculating shock price elasticity for the first type');
tic;
phi0 = ones(nx,1);

exp1 = CEimfd1(xx,tt,A,B,C,phi0,[X0;X0+dx;X0-dx],bc.SC);

% Shock elasiticy
tmp1 = exp1((nX0+1):(2*nX0),:);
tmp2 = exp1((2*nX0+1):(3*nX0),:);
exp1 = exp1(1:nX0,:);
tmpd = (tmp1-tmp2)/2/dx;
tmpld = tmpd ./ exp1; % #X0 by tt+1

nusigmax = diag( nux * sigmax.' ); %( #X0 by #shocks ) * (#X0 by #shock )'
nualphax = diag( nux * alphax ); % #X0 by 1
nusigmax = repmat(nusigmax,1,max(size(tt))+1); % #X0 by tt+1
nualphax = repmat(nualphax,1,max(size(tt))+1); % #X0 by tt+1

elas1.sg = nusigmax.*tmpld+nualphax;
toc;

%% Step 3.2 Second Type

disp('Calculating shock price elasticity for the second type');
tic;

phi0 = diag(Nu(xx) * alphaSC(xx));
exp2 = CEimfd1(xx,tt,A,B,C,phi0,X0,bc.SC);
elas2.sg = exp2 ./ exp1;
    
toc;

%% Step 4 Calculate final elasticity

elas1.out = elas1.g - elas1.sg;
elas2.out = elas2.g - elas2.sg;

end