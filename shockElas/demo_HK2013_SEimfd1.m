%% Step 0 Set up

clear all;
close all;
modelname = 'HK2013';
addpath( fullfile(modelname) );

%% Section 1 Model Construction

param.g = 0.02; % growth rate of endowment
param.sigma = 0.09; % volatility of endowment
param.m = 4; % tightness of financial friction (households willing to contribute maximum of m times expert wealth into the intermediary)
param.lambda = 0.6; % fraction of household assets in riskless bonds
param.rho = 0.04; % common discount rate
param.l = 1.84; % labor income to capital income ratio
param.xstar = (1-param.lambda) / (1-param.lambda+param.m); % point below which contraint starts binding (e.g., Xt <= xstar)
param.alphaX = @(x) 1./(1-param.lambda*(1-x)).*(x>param.xstar) + 1/(1+param.m)./x.*(x<=param.xstar); % expert leverage

model.muX = @(x) x.*(-param.l/(1+param.l).*param.rho+(param.alphaX(x)-1).^2.*param.sigma.^2);
model.sigmaX = @(x) x.*(param.alphaX(x)-1).*param.sigma;
model.muS = @(x) -param.rho/(1+param.l)-param.g+param.alphaX(x).*param.sigma.^2-0.5.*param.alphaX(x).^2*param.sigma.^2;
model.sigmaS = @(x) -param.alphaX(x).*param.sigma;

domain.x = [0 1];
domain.dt = 0.25;
domain.T = 100;
domain.nx = 10000;

%% Section 2 Get Density

% domainx = [0 1];
% BoundaryType(model,domainx,param.xstar,param);
[den,cumden] = Density(model,linspace(1e-5, 1-1e-5,10000),param.xstar,param,[],modelname);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(den.x,den.d,'linewidth',2); 
title('stationary density of $$X_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
if (~isempty(modelname))
    figname = strcat(modelname,'/diagnostic/density.eps');
    eval(['print -depsc2 ' figname])
end

figure; plot(cumden.x,cumden.cd); ylim([0 1]);
title('State CDF');

%% Section 3 X and SDF

tmpx = linspace(0,1,1000);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(tmpx,model.muX(tmpx),'--r','linewidth',2);
plot(tmpx,model.sigmaX(tmpx),'-b','linewidth',2);
title('dynamics of $$X_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
h=legend('$$\mu_X(x)$$','$$\sigma_X(x)$$','location','northeast');
set(h,'interpreter','latex');
ylim([-inf 0.1]);
if (~isempty(modelname))
    figname = strcat(modelname,'/diagnostic/muX_sigmaX.eps');
    eval(['print -depsc2 ' figname])
end

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(tmpx,-model.muS(tmpx)-0.5*model.sigmaS(tmpx).*model.sigmaS(tmpx),'--r','linewidth',2);
plot(tmpx,-model.sigmaS(tmpx),'-b','linewidth',2);
title('stochastic discount factor','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
h=legend('$$r(x)$$','$$\sigma_S(x)$$','location','northeast');
set(h,'interpreter','latex');
ylim([-0.2 0.8]);
if (~isempty(modelname))
    figname = strcat(modelname,'/diagnostic/muS_sigmaS.eps');
    eval(['print -depsc2 ' figname])
end

%% Section 4.1 Shock elasticities for first cash flow

model.muC = @(x) param.g - 0.5*param.sigma^2 + 0*x;
model.sigmaC = @(x)  param.sigma + 0*x; 

bc.C.l.const = @(t) 0;
bc.C.l.level = @(t) 0;
bc.C.l.deriv1 = @(t) 1;
bc.C.l.deriv2 = @(t) 0;
bc.C.u.const = @(t) 0;
bc.C.u.level = @(t) 0;
bc.C.u.deriv1 = @(t) 1;
bc.C.u.deriv2 = @(t) 0;
bc.SC.l.const = @(t) 0;
bc.SC.l.level = @(t) 0;
bc.SC.l.deriv1 = @(t) 1;
bc.SC.l.deriv2 = @(t) 0;
bc.SC.u.const = @(t) 0;
bc.SC.u.level = @(t) 0;
bc.SC.u.deriv1 = @(t) 1;
bc.SC.u.deriv2 = @(t) 0;

X0 = spline(cumden.cd(cumden.cd<0.99),cumden.x(cumden.cd<0.99),[0.1,0.5,0.9]).';
Nu = @(x) 1 + 0*x;
[tt, pr_el1_xx, pr_el2_xx] = SEimfd1(model,domain,bc,param,X0,Nu,[]);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el1_xx.g(1,:),'--r','linewidth',2);
plot([0 tt],pr_el1_xx.g(2,:),'-b','linewidth',2);
plot([0 tt],pr_el1_xx.g(3,:),':r','linewidth',2);
title('first-type shock-exposure elasticity for $$C_t$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','northeast');
ylim([0 0.15]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf1/elas1_exposure.eps');
    eval(['print -depsc2 ' figname])
end

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el1_xx.out(1,:),'--r','linewidth',2);
plot([0 tt],pr_el1_xx.out(2,:),'-b','linewidth',2);
plot([0 tt],pr_el1_xx.out(3,:),':r','linewidth',2);
title('first-type shock-price elasticity for $$C_t$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','southeast');
ylim([0 0.4]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf1/elas1.eps');
    eval(['print -depsc2 ' figname])
end

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el2_xx.g(1,:),'--r','linewidth',2);
plot([0 tt],pr_el2_xx.g(2,:),'-b','linewidth',2);
plot([0 tt],pr_el2_xx.g(3,:),':r','linewidth',2);
title('second-type shock-exposure elasticity for $$C_t$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','northeast');
ylim([0 0.15]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf1/elas2_exposure.eps');
    eval(['print -depsc2 ' figname])
end

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el2_xx.out(1,:),'--r','linewidth',2);
plot([0 tt],pr_el2_xx.out(2,:),'-b','linewidth',2);
plot([0 tt],pr_el2_xx.out(3,:),':r','linewidth',2);
title('second-type shock-price elasticity for $$C_t$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','southeast');
ylim([0 inf]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf1/elas2.eps');
    eval(['print -depsc2 ' figname])
end

%% Section 4.2 Shock elasticities for second cash flow

model.muC = @(x) -model.muS(x) - param.rho;
model.sigmaC = @(x) -model.sigmaS(x);

bc.C.l.const = @(t) 0;
bc.C.l.level = @(t) 0;
bc.C.l.deriv1 = @(t) 1;
bc.C.l.deriv2 = @(t) 0;
bc.C.u.const = @(t) 0;
bc.C.u.level = @(t) 0;
bc.C.u.deriv1 = @(t) 1;
bc.C.u.deriv2 = @(t) 0;
bc.SC.l.const = @(t) 0;
bc.SC.l.level = @(t) 0;
bc.SC.l.deriv1 = @(t) 1;
bc.SC.l.deriv2 = @(t) 0;
bc.SC.u.const = @(t) 0;
bc.SC.u.level = @(t) 0;
bc.SC.u.deriv1 = @(t) 1;
bc.SC.u.deriv2 = @(t) 0;

X0 = spline(cumden.cd(cumden.cd<0.99),cumden.x(cumden.cd<0.99),[0.1,0.5,0.9]).';
Nu = @(x) 1 + 0*x;
[tt, pr_el1_xx, pr_el2_xx] = SEimfd1(model,domain,bc,param,X0,Nu,[]);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el1_xx.g(1,:),'--r','linewidth',2);
plot([0 tt],pr_el1_xx.g(2,:),'-b','linewidth',2);
plot([0 tt],pr_el1_xx.g(3,:),':r','linewidth',2);
title('first-type shock-exposure elasticity for $$C_{e,t}$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','northeast');
ylim([0 inf]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf2/elas1_exposure.eps');
    eval(['print -depsc2 ' figname])
end

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el1_xx.out(1,:),'--r','linewidth',2);
plot([0 tt],pr_el1_xx.out(2,:),'-b','linewidth',2);
plot([0 tt],pr_el1_xx.out(3,:),':r','linewidth',2);
title('first-type shock-price elasticity for $$C_{e,t}$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','northeast');
ylim([0 inf]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf2/elas1.eps');
    eval(['print -depsc2 ' figname])
end

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el2_xx.g(1,:),'--r','linewidth',2);
plot([0 tt],pr_el2_xx.g(2,:),'-b','linewidth',2);
plot([0 tt],pr_el2_xx.g(3,:),':r','linewidth',2);
title('second-type shock-exposure elasticity for $$C_{e,t}$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','southeast');
ylim([0 inf]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf2/elas2_exposure.eps');
    eval(['print -depsc2 ' figname])
end

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot([0 tt],pr_el2_xx.out(1,:),'--r','linewidth',2);
plot([0 tt],pr_el2_xx.out(2,:),'-b','linewidth',2);
plot([0 tt],pr_el2_xx.out(3,:),':r','linewidth',2);
title('second-type shock-price elasticity for $$C_{e,t}$$','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','southeast');
ylim([0 inf]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf2/elas2.eps');
    eval(['print -depsc2 ' figname])
end

%% Section 5 Term Structure

bc.S.l.const = @(t) 0;
bc.S.l.level = @(t) 0;
bc.S.l.deriv1 = @(t) 1;
bc.S.l.deriv2 = @(t) 1;
bc.S.u.const = @(t) 0;
bc.S.u.level = @(t) 0;
bc.S.u.deriv1 = @(t) 1;
bc.S.u.deriv2 = @(t) 1;

% Notice that here I am using the boundary condition 0 = phi'(x) +
% phi''(x). This is because if just use 0 = phi'(x), then the solved term
% structure will be blowed up. I justify this by simulating a long term
% series of log S_t. The long-term converging point of the
% term-structure seems equal to the trend component of the log S_t
% process. Indeed, if one tries boundary condition 0 = phi(x), the results
% are exactly the same.

[tt, ts] = TSimfd1(model,domain,bc.S,param,X0,[]);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(tt,ts(1,:),'--r','linewidth',2);
plot(tt,ts(2,:),'-b','linewidth',2);
plot(tt,ts(3,:),':r','linewidth',2);
title('term structure of interest rates','interpreter','latex','fontsize',10);
xlabel('maturity (years)','interpreter','latex','fontsize',10);
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','southeast');
ylim([0 inf]);
if (~isempty(modelname))
    figname = strcat(modelname,'/diagnostic/ts.eps');
    eval(['print -depsc2 ' figname])
end

%% Section 6 Decomposition

process.muX = model.muX;
process.sigmaX = model.sigmaX;
process.muM = @(x)   model.muS(x);
process.sigmaM = @(x)   model.sigmaS(x);
domain.x = [0,1];
bc.M = bc.S;
[eig.M.eta,eig.M.e0,eig.M.xxe] = Decompose(process,domain,bc.M,param);
eig.M.e = @(x) spline(eig.M.xxe,eig.M.e0,x);   % col vec
eig.M.dloge = @(x) (   log( spline(eig.M.xxe,eig.M.e0,x+1e-4) ) - log( spline(eig.M.xxe,eig.M.e0,x-1e-4) )    ) / 2e-4;  % col vec

process.muX = @(x) model.muX(x) + model.sigmaX(x).*model.sigmaX(x).*eig.M.dloge(x) + model.sigmaX(x).*process.sigmaM(x);
denM = Density(process,linspace(0.001, 0.999,10000),param.xstar,param,[],[]);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(den.x,den.d,'linewidth',2);
plot(denM.x,denM.d,'--r','linewidth',2);
title('change-of-measure using $$S_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
if (~isempty(modelname))
    figname = strcat(modelname,'/diagnostic/density_changed_by_S.eps');
    eval(['print -depsc2 ' figname])
end

tmpx = linspace(0.001,1,1000);
gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(tmpx,model.muX(tmpx),'-b','linewidth',2);
plot(tmpx,process.muX(tmpx),'--r','linewidth',2);
title('change-of-measure using $$S_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
h=legend('$$\mu_X(x)$$','$$\tilde{\mu}_X(x)$$','location','northeast');
set(h,'interpreter','latex');
ylim([-inf 0.1]);
if (~isempty(modelname))
    figname = strcat(modelname,'/diagnostic/muX_changed.eps');
    eval(['print -depsc2 ' figname])
end

model.muC = @(x) param.g - 0.5*param.sigma^2 + 0*x;
model.sigmaC = @(x)  param.sigma + 0*x; 
process.muM = @(x) model.muC(x) + model.muS(x);
process.sigmaM = @(x) model.sigmaC(x) + model.sigmaS(x);
bc.M = bc.SC;
[eig.M.eta,eig.M.e0,eig.M.xxe] = Decompose(process,domain,bc.M,param);
eig.M.e = @(x) spline(eig.M.xxe,eig.M.e0,x);   % col vec
eig.M.dloge = @(x) (   log( spline(eig.M.xxe,eig.M.e0,x+1e-4) ) - log( spline(eig.M.xxe,eig.M.e0,x-1e-4) )    ) / 2e-4;  % col vec

process.muX = @(x) model.muX(x) + model.sigmaX(x).*model.sigmaX(x).*eig.M.dloge(x) + model.sigmaX(x).*process.sigmaM(x);
denM = Density(process,linspace(0.001, 0.999,10000),param.xstar,param,[],[]);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(den.x,den.d,'linewidth',2);
plot(denM.x,denM.d,'--r','linewidth',2);
title('change-of-measure using $$S_t C_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf1/density_changed.eps');
    eval(['print -depsc2 ' figname])
end

tmpx = linspace(0.001,1,1000);
gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(tmpx,model.muX(tmpx),'-b','linewidth',2);
plot(tmpx,process.muX(tmpx),'--r','linewidth',2);
title('change-of-measure using $$S_t C_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
h=legend('$$\mu_X(x)$$','$$\tilde{\mu}_X(x)$$','location','northeast');
set(h,'interpreter','latex');
ylim([-inf 0.1]);
if (~isempty(modelname))
    figname = strcat(modelname,'/cf1/muX_changed.eps');
    eval(['print -depsc2 ' figname])
end
