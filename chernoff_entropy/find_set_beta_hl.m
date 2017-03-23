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
hl = 120;
n_mu_y = 1;
%n_mu_z = 80;
n_mu_z = 1;
%n_beta = 80;
n_beta = 1;
%n_kappa = 60;
n_kappa = 600;
kappa_max = 0.16;
kappa_min = -0.16;
%xx = linspace(-0.05, mu_z_hat, n_mu_z);
%yy= linspace(beta_hat, 1.2, n_beta);
xx=linspace(mu_y_hat, mu_y_hat, n_mu_z);
yy = linspace(mu_z_hat, mu_z_hat, n_mu_z);
%yy= linspace(beta_hat, beta_hat, n_beta);
zz1 = linspace(kappa_min ,kappa_max,n_kappa);
zz2 = linspace(kappa_min ,kappa_max, n_kappa);

%t = chernoff_half_life_max(param, param_hat, hl);
for i = 1: n_mu_y
    tic;
    for j = 1:n_mu_z
        
        parfor k = 1:n_kappa
            [v1,FVAL,EXITFLAG]  = fminbnd( @(x) chernoff_half_life_max([xx(i) yy(j) x zz1(k)],param_hat, hl),...
                -beta_hat, beta_hat);
            if chernoff_half_life_max([xx(i) yy(j)  v1 zz1(k)],param_hat, hl) > .1
                v1 = 1;
            end
            v2 = fminbnd( @(x) chernoff_half_life_max([xx(i) yy(j) x zz1(k)],param_hat, hl),...
                beta_hat,2*beta_hat);
            if chernoff_half_life_max([xx(i) yy(j)  v2 zz1(k)],param_hat, hl) > .1
                v2 = 1;
            end
            v3 = fminbnd( @(x) chernoff_half_life_max([xx(i) yy(j) x zz1(k)],param_hat, hl),...
                -2.*beta_hat,  2.*beta_hat);
            if chernoff_half_life_max([xx(i) yy(j)  v3 zz1(k)],param_hat, hl) > .1
                v3 = 1;
            end
         
            vv = [v1  v2  v3 ];
            mm1 = max( vv);%, v2)
            vv(vv<0) = [];
            mm2 = min(vv );%, v2)
            if (mm2 == mm1)
                mm2 =10000;
            end
            ff1(i,j,k) = mm1;
            ff2(i,j,k) = mm2;
            disp([num2str(k)])
            disp( ff1(i,j,k) )
            disp( ff2(i,j,k) )
        end
        ff1 = squeeze(ff1);
        ff2 = squeeze(ff2);
        zz2(ff2==10000) = [];
        ff2(ff2==10000) = [];
        zz2(ff2==1) = [];
        ff2(ff2==1) = [];
         zz1(ff1==1) = [];
        ff1(ff1==1) = [];
        %         parfor k = 1:n_kappa
        %             v1 = fminbnd( @(x) chernoff_half_life_max([xx(i) yy(j) x zz1(k)],param_hat, hl),...
        %                 beta_hat+-0.0000001,1.6*beta_hat);
        %             mm = max(v1);%, v2);
        %             if chernoff_half_life_max([xx(i) yy(j)  mm zz1(k)],param_hat, hl) < .02
        %                 ff1(i,j,k) = mm;
        %             else
        %                 ff1(i,j,k) = 10000;
        %             end
        %             disp([num2str(k)])
        %             disp( ff1(i,j,k) )
        %         end
        %          zz1(ff1==10000) = [];
        %         ff1(ff1==10000) = [];
        disp(['i:' num2str(i) ' j:' num2str(j)]);
    end
    
    toc;
end

delete(par_obj)
save('beta_chernoff_set_120.mat')
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
