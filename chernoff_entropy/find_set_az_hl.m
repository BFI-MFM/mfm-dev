%clear all; close all;
clear all; close all;
c = parcluster('local');
c.NumWorkers = 16;
par_obj= parpool(c, c.NumWorkers);
if par_obj.NumWorkers>0

    display(['working in PARALLEL mode with ',num2str(matlabpool('size')),' workers']);

else
    profile clear
    profile on
    display('working in SERIAL mode');
end

tic;
mu_y_hat = 0.386;
mu_z_hat = 0;
beta_hat = 1;
kappa_hat = 0.019;
param_hat = [mu_y_hat  mu_z_hat beta_hat kappa_hat];

mu_z = mu_z_hat ;beta = beta_hat;
mu_y = mu_y_hat ;kappa = kappa_hat;

param = [mu_y  mu_z beta kappa];
hl = 360;
n_mu_y = 1;
%n_mu_z = 80;
n_mu_z = 1;
%n_beta = 80;
n_beta = 1;
%n_kappa = 60;
n_kappa = 500;
kappa_max = 0.06;
%xx = linspace(-0.05, mu_z_hat, n_mu_z);
%yy= linspace(beta_hat, 1.2, n_beta);
xx=linspace(mu_y_hat, mu_y_hat, n_mu_y);
%yy = linspace(mu_z_hat, mu_z_hat, n_mu_z);
yy= linspace(beta_hat, beta_hat, n_beta);
zz = linspace(0.000001,kappa_max, n_kappa);

 %t = chernoff_half_life_max(param, param_hat, hl);
 for i = 1: n_mu_y
     tic;
     for j = 1:n_beta
         parfor k = 1:n_kappa
            v1 = fminbnd( @(x) chernoff_half_life_max([xx(i) x  yy(j)  zz(k)],param_hat, hl),...
               -0.1, mu_z_hat);
             v2 = fminbnd( @(x) chernoff_half_life_max([xx(i) x yy(j)   zz(k)],param_hat, hl),...
               -0.1, mu_z_hat);
            ff(i,j,k) =min(v1, v2);
            disp([num2str(k)])
            disp( ff(i,j,k) )
         end
           disp(['i:' num2str(i) ' j:' num2str(j)]);
     end
     toc;
 end
 
delete(par_obj)
save('muz_chernoff_set_360.mat')
% 
% d = [zz(1)  zz(end)];
% nrank = 7;
% f= polyfit(zz, ff,nrank,domain(d));
% LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
% plot(zz, feval(f,zz)-ff,'r',LW,1.6);
% 
% dom_z = [-1 1];
% v = chebfun(@(z) z.^3, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
% z = chebfun(@(z) z, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
% u = feval(diff(f,1),zz)-diff(v,1).*z;
% rr = roots(u);
