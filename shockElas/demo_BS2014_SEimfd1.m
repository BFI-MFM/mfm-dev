%% Section 0 Set up

clear all;
close all;
modelname = 'BS2014';
addpath( fullfile(modelname) );

%% Section 1.1 Define parameters

param.a = 0.11; % productivity of experts
param.abar = 0.07; % productivity of households
param.delta = 0.05; % depreciation on capital
param.phi = 2; % adjustment cost parameter for cost function phi^(-1)*[exp(phi*investment_rate)-1]
param.sigma = 0.1; % fundamental shock volatility
param.rho = 0.07; % expert discount rate
param.rhobar = 0.04; % household discount rate (needs to be lower for stationarity)

%% Section 1.2 Solving for q(x)

tosolve = @(x,q) F(x,q,param);
param.q0 = (param.abar + 1/param.phi) ./ (param.rhobar + 1/param.phi);
param.dq0 = ( 1/param.sigma^2 * (param.a-param.abar) / param.q0 + 1 - (param.abar+1/param.phi)*(param.rho-param.rhobar)/(param.a-param.abar)/(param.rhobar+1/param.phi) ) * (param.a-param.abar) / (param.rhobar+1/param.phi);
sol = ode23s(tosolve, [1e-4 0.9], param.q0+1e-4*param.dq0);
plot(sol.x,sol.y);

param.Q2 = @(x) (param.a+1/param.phi) ./ ( (1-x)*param.rhobar+x*param.rho+1/param.phi ); % x >= param.xstar
hold on;
plot(sol.x,param.Q2(sol.x));

param.qq0 = [param.q0 sol.y(sol.x<0.5)];
param.xx0 = [0 sol.x(sol.x<0.5)];
plot(param.xx0,param.qq0);

param.Q1 = @(x) spline(param.xx0,param.qq0,x); % x < param.xstar
param.xstar = fzero(@(x) param.Q1(x)-param.Q2(x),0.3);
plot(param.xstar,param.Q1(param.xstar),'*');

%% Section 1.3 dQ polish

param.eps = 1e-4;

param.dQ1 = @(x) ( ( param.Q1(x+param.eps) - param.Q1(x-param.eps) ) / 2/param.eps ) ;
xx = linspace(0.001,0.5,100);
figure;plot(xx,param.dQ1(xx));

xx1 = linspace(0.001,0.8,1000); % 0.35 for previous plots
qq1 = param.dQ1(xx1);
param.xx1 = [0 xx1];
param.qq1 = [param.dq0 qq1];
param.dQ1 = @(x) spline(param.xx1,param.qq1,x) ;
xx = linspace(0,0.4,100);
hold on; plot(xx,param.dQ1(xx));

param.d2Q1 = @(x) ( ( param.dQ1(x+param.eps) - param.dQ1(x-param.eps) ) / 2/param.eps ) ;
xx = linspace(0,0.5,100);
figure;plot(xx,param.d2Q1(xx));

param.dQ2 = @(x) ( param.Q2(x+param.eps) - param.Q2(x-param.eps) ) / 2/param.eps; % x >= param.xstar
param.d2Q2 = @(x) ( param.Q2(x+param.eps) + param.Q2(x-param.eps) - 2*param.Q2(x)) / param.eps^2; % x >= param.xstar

%% Section 1.4 Model Construction

param.Q = @(x) param.Q1(x) .* ( x<param.xstar) + param.Q2(x) .* (x>param.xstar);
param.dQ = @(x) param.dQ1(x) .* ( x<param.xstar) + param.dQ2(x) .* (x>param.xstar);
param.d2Q = @(x) param.d2Q1(x) .* ( x<param.xstar) + param.d2Q2(x) .* (x>param.xstar);

model.sigmaX = @(x) sigmaXbs(x,param);
model.muX = @(x) muXbs(x,param);
model.muS = @(x) - ( ...
    param.rho + 1/param.phi*log(param.Q(x)) - param.delta ...
    + param.dQ(x)./param.Q(x).*model.muX(x) + 0.5*param.d2Q(x)./param.Q(x).*model.sigmaX(x).^2 ...
    + model.muX(x)./x - 0.5 * ( param.sigma^2 + ( param.dQ(x)./param.Q(x).*model.sigmaX(x) ).^2 + (model.sigmaX(x)./x).^2 ) ...
    );
model.sigmaS = @(x) - ( ...
    param.sigma + ( param.dQ(x)./param.Q(x).*model.sigmaX(x) ) + model.sigmaX(x)./x ...
    );

domain.x = [0 1];
domain.dt = 0.25;
domain.T = 100;
domain.nx = 20000;

%% Section 2 Get Density

% domainx = [0.00000001 0.99999999];
% BoundaryType(model,domainx,param.xstar,param);
[den,cumden] = Density(model,linspace(1e-5, 1-1e-5,10000),param.xstar,param,0.999,modelname);

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
plot(den.x,den.d,'linewidth',2); 
title('stationary density of $$X_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
ylim([0 3]);
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

gcf = figure;
set(gcf, 'Units', 'Inches', 'Position', [0,0,5,3], 'PaperUnits', 'Inches', 'PaperPosition', [0,0,5,3])
hold on; box on;
tmp = param.dQ(tmpx)./param.Q(tmpx).*model.muX(tmpx) + 0.5*param.d2Q(tmpx)./param.Q(tmpx).*model.sigmaX(tmpx).^2;
plot(tmpx,tmp,'--r','linewidth',2);
tmp = param.dQ(tmpx)./param.Q(tmp).*model.sigmaX(tmpx) ;
plot(tmpx,tmp,'-b','linewidth',2);
title('dynamics of $$q_t$$','interpreter','latex','fontsize',10);
xlabel('expert wealth share $$x$$','interpreter','latex','fontsize',10);
h=legend('$$\mu_q(x)$$','$$\sigma_q(x)$$','location','northeast');
set(h,'interpreter','latex');
ylim([-inf inf]);
if (~isempty(modelname))
    figname = strcat(modelname,'/diagnostic/muq_sigmaq.eps');
    eval(['print -depsc2 ' figname])
end

%% Section 4.1 Shock elasticities for first cash flow

model.muC = @(x) 1./param.phi.*log(param.Q(x)) - param.delta ...
    + param.dQ(x)./param.Q(x).*model.muX(x) + 0.5*param.d2Q(x)./param.Q(x).*model.sigmaX(x).^2 ...
    + (param.rho-param.rhobar)./((1-x).*param.rhobar+x.*param.rho).*model.muX(x) ...
    - 0.5 * ( param.sigma^2 + ( param.dQ(x)./param.Q(x).*model.sigmaX(x) ).^2 + ( (param.rho-param.rhobar)./((1-x).*param.rhobar+x.*param.rho) ).^2 .* (model.sigmaX(x)).^2 ) ;
model.sigmaC = @(x) param.sigma + ( param.dQ(x)./param.Q(x).*model.sigmaX(x) ) + model.sigmaX(x).*(param.rho-param.rhobar)./((1-x).*param.rhobar+x.*param.rho);

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
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','east');
ylim([0 inf]);
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
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','east');
ylim([0 inf]);
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
legend('x = 10th percentile','x = 50th percentile','x = 90th percentile','location','southeast');
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

model.muC = @(x) - model.muS(x) - param.rho;
model.sigmaC = @(x) - model.sigmaS(x);

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
bc.S.l.deriv2 = @(t) 0;
bc.S.u.const = @(t) 0;
bc.S.u.level = @(t) 0;
bc.S.u.deriv1 = @(t) 1;
bc.S.u.deriv2 = @(t) 0;

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
process.muM = @(x)  model.muS(x);
process.sigmaM = @(x) model.sigmaS(x);
[eig.S.eta,eig.S.e0,eig.S.xxe] = Decompose(process,domain,bc.S,param);
eig.S.e = @(x) spline(eig.S.xxe,eig.S.e0,x);   % col vec
eig.S.dloge = @(x) (   log( spline(eig.S.xxe,eig.S.e0,x+1e-4) ) - log( spline(eig.S.xxe,eig.S.e0,x-1e-4) )    ) / 2e-4;  % col vec

process.muX = @(x) model.muX(x) + model.sigmaX(x).*model.sigmaX(x).*eig.S.dloge(x) + model.sigmaX(x).*process.sigmaM(x);
denS = Density(process,linspace(0.0001, 0.9999,10000),param.xstar,param,[],[]);

figure;plot(den.x,den.d);
hold on;plot(denS.x,denS.d); ylim([0,max(den.d)]);
title('Changed Density');
saveas(gcf,strcat(modelname,'/diagnostic/density_changed_by_S.png'));

figure;plot(tmpx,process.muX(tmpx));hold on;
plot(tmpx,process.sigmaX(tmpx));
title('$\mu_X (x)$ and $\sigma_X (x)$','Interpreter','latex');
saveas(gcf,strcat(modelname,'/diagnostic/','muX_sigmaX_changed.png'));
