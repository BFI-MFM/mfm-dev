tic;
%setup parameters (latest specs)
param.r = 0.0;
param.rho = 0.05;
param.sigma0 = 0.1;
param.p = 0.02;
param.beta = 1;
param.R = 0.3;
param.gamma = 0.5;

LW = 'linewidth'; lw = 1.6;  FS = 'fontsize';
%arbitrary starting value for the integral
int_R = 0;
%error toletance
tol_E = 1e-03;
%arbitrary starting value for upper bound of E
upper_E = 0.1;
mc_equity = log(1+param.gamma);
alpha = (param.R-param.p).^(-param.beta);
Lfun = @(r) alpha*(param.R-r).^param.beta;
diff_Lfun = @(r) -alpha.*param.beta.*(param.R-r).^(param.beta-1);
while abs(int_R-mc_equity) > tol_E
    domain_E = [0 upper_E];
    Rfun = rochet_fun(domain_E, param);
    int_R = sum((Rfun-param.p-param.r)./param.sigma0^2./Lfun(Rfun));
    error_E = int_R-mc_equity;
    %update step
    upper_E = upper_E - 0.1*error_E./mc_equity;
    %disp(int_R);  disp(upper_E)
end
Rfun = rochet_fun(domain_E, param);
diff_Rfun = @(e,r)	-1/param.sigma0^2.*(2*(param.rho-param.r).*param.sigma0^2+ ...
    (r-param.p-param.r).^2+2*param.r*e.*(r-param.p-param.r)./Lfun(r))./...
    (Lfun(r)-diff_Lfun(r).*(r-param.p-param.r));
xfun = chebfun(@(x) x,domain_E);

u = @(e) sum((Rfun-param.p-param.r)./param.sigma0^2./Lfun(Rfun), e, domain_E(2));

sigmacheb = -param.sigma0.*Lfun(Rfun);
mucheb = param.r*xfun+Lfun(Rfun).*(Rfun-param.p-param.r);
pull = mucheb./sigmacheb-diff(sigmacheb)/2.;

xi = diff_Rfun(xfun,Rfun)./Rfun;
betacheb = mucheb.*xi+1/2*sigmacheb.^2.*diff(xi,1);
alphacheb = sigmacheb.*xi;

mbook = chebfun(u, domain_E,'splitting','on', 'vectorize');
xi_sdf = diff(mbook,1);
betacheb_sdf = -param.rho+ mucheb.*xi_sdf+1/2*sigmacheb.^2.*diff(xi_sdf,1);
alphacheb_sdf = sigmacheb.*xi_sdf;

% plot on grids
ne = 100;
ee = linspace(domain_E(1), domain_E(end), ne);
for i = 1:ne
    uu(i) = feval(u,ee(i));
end

figure;
subplot(2,1,1); hold on
plot(ee,feval(mucheb,ee),LW, lw);
xlabel('E', FS,16)
ylabel('\mu(E)',FS,16);

subplot(2,1,2); hold on
plot(ee, feval(sigmacheb,ee),LW, lw);
xlabel('E',FS,16)
ylabel('\sigma(E)',FS,16);

figure;
subplot(2,1,1); hold on
plot(ee, exp(uu),LW, lw);
xlabel('E', FS,16)
ylabel('market-to-book',FS,16);

subplot(2,1,2); hold on
plot(ee, feval(Rfun,ee),LW, lw);
xlabel('E',FS,16)
ylabel('R(E)',FS,16);
domain_x = domain_E;
l.boundary_type ='reflecting';
r.boundary_type ='reflecting';
m_name ='interest rate';
sh_type ='state space: bank equity';
toc;