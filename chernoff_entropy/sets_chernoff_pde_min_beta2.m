clear all; close all;
% load beta_chernoff_set_960.mat;
% ff = squeeze(ff)';
% zz(ff==min(ff))=[];
% ff(ff==min(ff))=[];
% ff1 =ff;
% zz1 =zz;
% load beta_chernoff_set_960.mat;
% ff = squeeze(ff)';
% zz(ff==min(ff))=[];
% ff(ff==min(ff))=[];
% ff2 =ff;
% zz2 =zz;
% load beta_chernoff_set_960.mat;
% ff = squeeze(ff)';
% zz(ff==min(ff))=[];
% ff(ff==min(ff))=[];
% ff3 =ff;
% zz3 =zz;
%  figure; plot(zz1,ff1,LW,1.6); hold on;
%  plot(zz2,ff2,LW,1.6); plot(zz3,ff3,LW,1.6);
% xlabel('\kappa', FS, 14); ylabel('\beta',FS, 14);
%  leg=legend('T=960', 'T=960', 'T=960');
%   set(leg,FS, 14); set(leg,'Interpreter','latex');
% print(gcf, '-depsc2', '-loose','chernoff_set_beta.eps')
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
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
% load muz_chernoff_set_960.mat;
% ff = squeeze(ff)';
% zz(ff==max(ff))=[];
% ff(ff==max(ff))=[];
% ff1 =ff;
% zz1 =zz;
% load muz_chernoff_set_960.mat;
% ff = squeeze(ff)';
% zz(ff==max(ff))=[];
% ff(ff==max(ff))=[];
% ff2 =ff;
% zz2 =zz;
%  load beta_chernoff_set_960.mat;
%  ff = squeeze(ff)';
%  zz(ff==min(ff))=[];
%  ff(ff==min(ff))=[];
%   figure;
% plot(zz,ff,LW,1.6);
%  xlabel('\kappa', FS, 14); ylabel('\beta',FS, 14);
%  print(gcf, '-depsc2', '-loose','chernoff_set_beta_960.eps')
load beta_chernoff_set_960.mat;
%ff1 = squeeze(ff1)';
%zz1(ff1==min(ff1))=[];
%ff1(ff1==min(ff1))=[];
%ff2 = squeeze(ff2)';
%zz2(ff2==max(ff2))=[];
%ff2(ff2==max(ff2))=[];
mu_y_hat = 0.386;
mu_z_hat =  0;
beta_hat = 1;
kappa_hat = 0.019;
sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';
delta = 0.002;
sigma = [sigma_y'; sigma_z'];
sigma_inv = inv(sigma);
ns = 60;
chernoff_s = linspace(0.1,0.9,ns);
d1 = [zz1(1)  zz1(end)];
nrank = 24;
f1= polyfit(zz1, ff1',nrank,domain(d1));
figure;plot(zz1, feval(f1,zz1)-ff1','r',LW,1.6);
figure; plot(f1,LW,1.6);
xlabel('\kappa', FS, 14); ylabel('\beta',FS, 14);
%print(gcf, '-depsc2', '-loose','chernoff_set_beta_960.eps')
d2 = [zz2(1)  zz2(end)];
f2= polyfit(zz2, ff2',nrank,domain(d2));
figure;plot(zz2, feval(f2,zz2)-ff2','r',LW,1.6);
figure; plot(f2,LW,1.6);
xlabel('\kappa', FS, 14); ylabel('\beta',FS, 14);
n_theta = 1;
%theta_space = linspace(0.08,2,n_theta);
%theta_space = linspace(0.08,2,n_theta);
theta_space = 0.1178;
%theta_space = 0.0937;
%theta_space = 10;
ss_sigma  = sigma*sigma';
dom_z = [-2 0 2];
xx = linspace(dom_z(1),dom_z(3),12);
onefun = chebfun(@(z) 1, dom_z, 'eps', 1e-6, 'vectorize','splitting', 'on');
zfun = chebfun(@(z) z, dom_z, 'eps', 1e-6, 'vectorize','splitting', 'on');
for kk=1:n_theta
    disp('computing theta');
    %kk=1;
    theta = theta_space(kk);
    disp(theta);
    %valfun = 0.05.*chebfun(@(z) z+cos(0.3.*z-1), dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
    valfun = chebfun(@(z) 0.01.*(0.3.*z-2).^3, dom_z, 'eps', 1e-6, 'vectorize');
    %valfun= w{3000};
    w{1} = valfun;
    tol_val = 1; i=1;
    while (tol_val > 0.01)
        i=i+1;
        %FOC conditon
        u1 = 0.01.*feval(diff(f1,1),zz1)-diff(valfun,1);
        temp = roots(u1);
        z_roots1 =temp(1,:);
        zzr1 = zz1;ffr1 = ff1';
        %if isnan(z_roots)
        zzr1(isnan(z_roots1)) = [];
        ffr1(isnan(z_roots1)) = [];
        z_roots1(isnan(z_roots1)) = [];
        %end
        dr1 = [min(zzr1) max(zzr1)];
        roots_fun1 = polyfit(zzr1,z_roots1,4,domain(dr1));
        %figure;plot(roots_fun);hold on; plot(zz2,z_roots);
        kappa_fun1 = inv(roots_fun1,'eps', 1e-5);
        aa = min(roots_fun1);
        if aa < dom_z(1)
            aa = dom_z(1)+0.0001;
        end
        u2 = 0.01.*feval(diff(f2,1),zz2)-diff(valfun,1);
        temp = roots(u2);
        z_roots2 =temp(1,:);
        zzr2 = zz2;ffr2 = ff2';
        %if isnan(z_roots)
        zzr2(isnan(z_roots2)) = [];
        ffr2(isnan(z_roots2)) = [];
        z_roots2(isnan(z_roots2)) = [];
        %end
        dr2 = [min(zzr2) max(zzr2)];
        roots_fun2 = polyfit(zzr2,z_roots2,4,domain(dr2));
        %figure;plot(roots_fun);hold on; plot(zz2,z_roots);
        kappa_fun2 = inv(roots_fun2,'eps', 1e-5);
        %if min(roots_fun) < dom_z(1)
        %	min(roots_fun) = dom_z(1);
        %end
        bb = max(roots_fun2);
        if bb > dom_z(end)
            bb = dom_z(end)-0.0001;
        end
        z_sol = chebfun({@(z) feval(kappa_fun1,aa),@(z) feval(kappa_fun1,z),@(z) feval(kappa_fun2,z), @(z) feval(kappa_fun2,bb) }, [dom_z(1) aa  0 bb dom_z(end)], ...
            'eps', 1e-4, 'vectorize');
        z_sol1 = chebfun({@(z) feval(kappa_fun1,aa),@(z) feval(kappa_fun1,z)}, [dom_z(1) aa  0], ...
            'eps', 1e-4, 'vectorize');
        z_sol2 = chebfun({@(z) feval(kappa_fun2,z), @(z) feval(kappa_fun2,bb)}, [0 bb dom_z(end)], ...
            'eps', 1e-4, 'vectorize');
        f1.domain = [min(z_sol1) max(z_sol1)];
        f2.domain =  [min(z_sol2) max(z_sol2)];
        f = chebfun({@(z) feval(f1,feval(z_sol1,z)), @(z) feval(f2,feval(z_sol2,z))} , [dom_z(1)  0  dom_z(end)], ...
            'eps', 1e-5, 'vectorize');
        %f.domain = z_sol.domain;
        N = chebop(dom_z);
        % N.op = 	@(v) -delta.*v + 0.01.*(f(z_sol)+beta.*zfun)+diff(v,1).*(mu_z-z_sol.*zfun)+...
        % 	1/2.*norm(sigma_z).^2.*diff(v,2) - theta_term(v);
        theta_term = @(v)  (0.01.*(0.01*ss_sigma(1,1)+diff(v,1).*ss_sigma(1,2))+...
            diff(v,1).*(0.01.*ss_sigma(2,1)+diff(v,1).*ss_sigma(2,2)))...
            *1./(2.*theta);
        N.op = 	@(v) -delta.*v + 0.01.*(mu_y+f.*zfun)+diff(v,1).*(mu_z-z_sol.*zfun)+...
            1/2.*norm(sigma_z).^2.*diff(v,2)- theta_term(v);
        delta_f = N(w{i-1});
        v =   w{i-1} + delta_f;
        %w{i} = v;
        w{i} = chebfun.spline(xx,v(xx));
        valfun = w{i};
        tol_val = (abs(max(delta_f ))+abs(min(delta_f )))./abs(max(v));
        disp(i); disp(tol_val);
    end
    disp(i); disp((abs(max(delta_f ))+abs(min(delta_f )))./abs(max(v)));
    %[feval(diff(v,1), -1) 0.01/(delta+min(z_sol))]
    %[feval(diff(v,1), 1) 0.01/(delta+max(z_sol))]
    %[feval(diff(v,1), 20) 0.01/(delta+max(z_sol))]
    figure; plot(f.*zfun.*0.01+ mu_y_hat.*0.01.*onefun, LW,1.6);
    xlabel('z', FS, 14); ylabel('0.01*\beta z+\mu_y',FS, 14);
    %print(gcf, '-depsc2', '-loose','sets_drift1.eps')
    figure; plot(diff(valfun,1), LW,1.6);
    xlabel('z', FS, 14); ylabel('d v/ dz',FS, 14);
    %print(gcf, '-depsc2', '-loose','sets_drift_val_derivative.eps')
    figure; plot(-z_sol.*zfun, LW,1.6);hold on;
    plot(-kappa_hat.*zfun, LW,1.6);
    xlabel('z', FS, 14);
    leg=legend('$-\kappa(z) z$','$-\hat{\kappa} z$', 'Location', 'best');
    set(leg,FS, 14); set(leg,'Interpreter','latex');
    %print(gcf, '-depsc2', '-loose','sets_beta_960_kappa_drift.eps')
    figure; plot(f, LW,1.6);hold on;
    plot(beta_hat.*onefun, LW,1.6);
    xlabel('z', FS, 14); leg=legend('$\beta(z)$','$\hat{\beta}$', 'Location', 'best');
    set(leg,FS, 14); set(leg,'Interpreter','latex');
    %print(gcf, '-depsc2', '-loose','beta_960_beta(z).eps')
    rfun1 = sigma_inv(1,1).*(mu_y-mu_y_hat)+sigma_inv(1,2).*(mu_z-mu_z_hat)+...
        sigma_inv(1,1).*(f-beta_hat).*zfun+sigma_inv(1,2).*(kappa_hat-z_sol).*zfun;
    rfun2 = sigma_inv(2,1).*(mu_y-mu_y_hat)+sigma_inv(2,2).*(mu_z-mu_z_hat)+...
        sigma_inv(2,1).*(f-beta_hat).*zfun+sigma_inv(2,2).*(kappa_hat-z_sol).*zfun;
    hfun1 =   rfun1 - 1./theta.*(sigma(1,1).*0.01+sigma(2,1).*diff(valfun,1));
    hfun2 =  rfun2 - 1./theta.*(sigma(1,2).*0.01+sigma(2,2).*diff(valfun,1));
    figure;plot(-hfun1,LW,1.6);hold on;plot(-hfun2,'r',LW,1.6); xlabel('z', FS, 14); ylabel('-h^*(z)',FS, 14);
    h=legend('-h_1(z)','-h_2(z)','Location', 'best'); set(h,FS, 14);
    %file_save = ['sets_density-theta-xi-0-01-' num2str(kk) '.mat'];
    %save(file_save);
    parfor jj = 1:ns
        N_c = chebop(dom_z);
        s = chernoff_s(jj);
        tic;
        N_c.op = 	@(e) -0.5.*s.*(1-s).*e.*(hfun1.^2+hfun2.^2) +...
            (s.*(sigma_z(1).*hfun1+sigma_z(2).*hfun2)+ mu_z_hat-kappa_hat.*zfun).*diff(e,1)+...
            1/2.*norm(sigma_z).^2.*diff(e,2);
        N_c.bc = 'neumann';
        [V,D]= eigs(N_c,100);
        e_v =  -sort((real(diag(D))),'descend');
        disp(e_v(1));
        rho1(jj,kk) = e_v(1);
        disp(jj);
        toc;
        N2_c = chebop(dom_z);
        s = chernoff_s(jj);
        tic;
        N2_c.op = 	@(e) -0.5.*s.*(1-s).*e.*(rfun1.^2+rfun2.^2) +...
            (s.*(sigma_z(1).*rfun1+sigma_z(2).*rfun2) + mu_z_hat-kappa_hat.*zfun).*diff(e,1)+...
            1/2.*norm(sigma_z).^2.*diff(e,2);
        N2_c.bc = 'neumann';
        [V,D]= eigs(N2_c,100);
        e_v =  -sort((real(diag(D))),'descend');
        disp(e_v(1));
        rho2(jj,kk) = e_v(1);
        disp(jj);
        toc;
    end
    file_save = ['hl_theta_beta_960_5-' num2str(kk) '.mat'];
    %file_save = ['hl_theta_beta_960_120.mat'];
    save(file_save);
end
delete(par_obj)