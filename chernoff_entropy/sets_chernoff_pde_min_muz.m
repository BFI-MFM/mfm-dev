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

LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';

set_path_tensor

  load muz_chernoff_set_960.mat;
  ff = squeeze(ff)';
  zz(ff==max(ff))=[];
  ff(ff==max(ff))=[];

    figure; plot(zz,ff,LW,1.6); hold on; 
% % %  plot(zz2,ff2,LW,1.6); 
% % plot(zz3,ff3,LW,1.6);
   xlabel('\kappa', FS, 14); ylabel('\alpha_z',FS, 14); 
% % %  leg=legend('T=120', 'T=960', 'T=360');
% % %   set(leg,FS, 14); set(leg,'Interpreter','latex');
   print(gcf, '-depsc2', '-loose','chernoff_set_alpha_z_960.eps')


load muz_chernoff_set_960.mat;

ff = squeeze(ff)';
zz(ff==max(ff))=[];
ff(ff==max(ff))=[];

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

d = [zz(1)  zz(end)];
nrank = 18;
f= polyfit(zz, ff,nrank,domain(d));
figure;plot(zz, feval(f,zz)-ff,'r',LW,1.6);

n_theta = 1;
%theta_space = [linspace(.15,.3,15) linspace(1,20,5)];
%theta_space = 0.3788;
%theta_space = 0.4563;
%theta_space = 0.4;
%theta_space = 0.3412;
%theta_space = 0.1932;
%theta_space = 0.2215;
theta_space = 0.2504;
ss_sigma  = sigma*sigma';
dom_z = [-1 0 1];
onefun = chebfun(@(z) 1, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
zfun = chebfun(@(z) z, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');

%FOC condition
u = feval(diff(f,1),zz)-zfun;

z_roots = roots(u);
zz2 = zz;
zz2(isnan(z_roots)) = [];
ff2 = ff;
ff2(isnan(z_roots)) = [];
z_roots(isnan(z_roots)) = [];
d2 = [min(zz2) max(zz2)];
roots_fun = polyfit(zz2,z_roots,nrank,domain(d2));
figure; plot(roots_fun);hold on; plot(zz2,z_roots);
kappa_fun = inv(roots_fun);

z_sol = chebfun({@(z) min(zz2),@(z) feval(kappa_fun,z), @(z) max(zz2) }, [dom_z(1) min(roots_fun) max(roots_fun) dom_z(end)], ...
			'eps', 1e-6, 'vectorize');

f.domain = z_sol.domain;

zzz= linspace(-1,1,100);
figure; plot(zzz,feval(f(z_sol),zzz)-feval(zfun.*z_sol,zzz),LW,1.6);hold on; plot(zzz,mu_z-feval(zfun.*kappa_hat,zzz),LW,1.6); 
xlabel('z'); leg=legend('$\alpha_z-z*\kappa$', '$\hat{\alpha_z}-z*\hat{\kappa}$'); 
 set(leg,FS, 14); set(leg,'Interpreter','latex');
%print(gcf, '-depsc2', '-loose','sets_implied_drift_alpha_z_960.eps')

for kk=1:n_theta
	disp('computing theta');
		%kk=1;
	theta = theta_space(kk);
	disp(theta);
	valfun = chebfun(@(z) z+cos(0.3.*z-1), dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
	%valfun= w{3000};
	zfun = chebfun(@(z) z, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
	w{1} = valfun;
	for i = 2:1000
        f.domain = z_sol.domain;
		N = chebop(dom_z);

 		 theta_term = @(v)  (0.01.*(0.01.*ss_sigma(1,1)+diff(v,1).*ss_sigma(1,2))+...
			 diff(v,1).*(0.01.*ss_sigma(2,1)+diff(v,1).*ss_sigma(2,2)))...
			 *1./(2.*theta);
		N.op = 	@(v) -delta.*v + 0.01.*(mu_y+beta.*zfun)+diff(v,1).*(f(z_sol)-z_sol.*zfun)+...
			1/2.*norm(sigma_z).^2.*diff(v,2)- theta_term(v);

		delta_f = N(w{i-1});
		v =   w{i-1} + delta_f;
		xx = linspace(dom_z(1),dom_z(3),12);
		%w{i} = v;
		w{i} = chebfun.spline(xx,v(xx));
		valfun = w{i};
		disp(i); disp((abs(max(delta_f ))+abs(min(delta_f )))./abs(max(v)));
	end
	disp(i); disp((abs(max(delta_f ))+abs(min(delta_f )))./abs(max(v)));

	%[feval(diff(v,1), -1) 0.01/(delta+min(z_sol))]
	%[feval(diff(v,1), 1) 0.01/(delta+max(z_sol))]
	%[feval(diff(v,1), 20) 0.01/(delta+max(z_sol))]

	% figure; plot(f(z_sol).*0.01-z_sol.*zfun, LW,1.6); 
	% xlabel('z', FS, 14); ylabel('0.01*\alpha_y-z\kappa(z)',FS, 14); 
	% print(gcf, '-depsc2', '-loose','sets_drift1.eps')

	% figure; plot(f(z_sol)+zfun, LW,1.6); 
	% xlabel('z', FS, 14); ylabel('\alpha_y+z',FS, 14); 
	% print(gcf, '-depsc2', '-loose','sets_drift2.eps')

	% 	figure; plot(diff(valfun,1), LW,1.6); 
	% xlabel('z', FS, 14); ylabel('d v/ dz',FS, 14); 
	% print(gcf, '-depsc2', '-loose','sets_drift_val_derivative.eps')

	% figure; plot(-z_sol.*zfun, LW,1.6);hold on;
	%  plot(-kappa_hat.*zfun, LW,1.6);
	%   xlabel('z', FS, 14); 
	% leg=legend('$-\kappa(z) z$','$-\hat{\kappa} z$', 'Location', 'best'); 
	% set(leg,FS, 14); set(leg,'Interpreter','latex');
	% print(gcf, '-depsc2', '-loose','sets_drift_hat.eps')

	% 	figure; plot(f(z_sol), LW,1.6);hold on;
	%  plot(mu_y_hat.*onefun, LW,1.6);
	%   xlabel('z', FS, 14); leg=legend('$\alpha_y(z)$','$\hat{\alpha_y}$', 'Location', 'best'); 
	%   set(leg,FS, 14); set(leg,'Interpreter','latex');
	% print(gcf, '-depsc2', '-loose','sets_drift_hat2.eps')

	rfun1 = sigma_inv(1,1).*(mu_y-mu_y_hat)+sigma_inv(1,2).*(f(z_sol)-mu_z_hat)+...
		sigma_inv(1,1).*(beta-beta_hat).*zfun+sigma_inv(1,2).*(kappa_hat-z_sol).*zfun;
	
	rfun2 = sigma_inv(2,1).*(mu_y-mu_y_hat)+sigma_inv(2,2).*(f(z_sol)-mu_z_hat)+...
		sigma_inv(2,1).*(beta-beta_hat).*zfun+sigma_inv(2,2).*(kappa_hat-z_sol).*zfun;
   
   	hfun1 =   rfun1 - 1./theta.*(sigma(1,1).*0.01+sigma(2,1).*diff(valfun,1));
	hfun2 =  rfun2 - 1./theta.*(sigma(1,2).*0.01+sigma(2,2).*diff(valfun,1));

    figure;plot(-hfun1,LW,1.6);hold on;plot(-hfun2,'r',LW,1.6); xlabel('z', FS, 14); ylabel('-h^*(z)',FS, 14); 
    h=legend('-h_1(z)','-h_2(z)','Location', 'best'); set(h,FS, 14); 
    %print(gcf, '-depsc2', '-loose','sets_h_alpha_z_960.eps')

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
 	%file_save = ['hl_theta_alpha_z_960_5-' num2str(kk) '.mat'];
 	file_save = ['hl_theta_alpha_z_960_120.mat'];
	save(file_save);

end 

delete(par_obj)


% load sets_h_960_T_80_muz.mat
% hfun1_80 = hfun1;
% hfun2_80 = hfun2;
% rfun1_80 = rfun1;
% rfun2_80 = rfun2;
% load sets_h_960_T_100_muz.mat
% hfun1_100 = hfun1;
% hfun2_100 = hfun2;
% rfun1_100 = rfun1;
% rfun2_100 = rfun2;

% zzz =linspace(-0.5,0.5, 500);

% FS = 'fontsize';
% set(gcf,'paperpositionmode','auto')
% f1= figure;
% plot(zzz,sigma_z(1).*feval(hfun1_80,zzz)+sigma_z(2).*feval(hfun2_80,zzz),'Linewidth',1.6');hold on;
% plot(zzz,sigma_z(1).*feval(hfun1_100,zzz)+sigma_z(2).*feval(hfun2_100,zzz),'Linewidth',1.6');
% plot(zzz,sigma_z(1).*feval(rfun1_80,zzz)+sigma_z(2).*feval(rfun2_80,zzz),'Linewidth',1.6');hold on;
% plot(zzz,sigma_z(1).*feval(rfun1_100,zzz)+sigma_z(2).*feval(rfun2_100,zzz),'Linewidth',1.6');
% xlabel('z',FS,14);
% leg = legend('$\sigma_z\cdot H, T=80$','$\sigma_z\cdot H, T=100 $',...
% 		'$\sigma_z\cdot R, T=80$','$\sigma_z\cdot R, T=100 $');
% set(leg,'Interpreter','latex');
% set(leg,FS,14, 'Location','best');
% %D = f1.PaperPosition; % Returns 1x4 vector [left right width height]
% %f1.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
% % print -f1 -dpsc2 sets_drift1.eps;
% print(gcf,'-depsc2','-loose','sets_alpha_z_960_h.eps');
% %print('-depsc2', 'sets_drift1.eps')


% FS = 'fontsize';
% set(gcf,'paperpositionmode','auto')
% figure;
% plot(zzz,feval(hfun1_80,zzz),'Linewidth',1.6');hold on;plot(zzz,feval(hfun1_100,zzz),'Linewidth',1.6');
% plot(zzz,feval(rfun1_80,zzz),'Linewidth',1.6');hold on;plot(zzz,feval(rfun1_100,zzz),'Linewidth',1.6');
% xlabel('z',FS,14);
% leg = legend('$h1, T=80$','$h1, T=100 $',...
% 		'$r1, T=80$','$r1, T=100 $' );
% set(leg,'Interpreter','latex');
% set(leg,FS,14, 'Location','best');
% print(gcf,'-depsc2','-loose','sets_alpha_z_960_h1.eps');
% figure;
% plot(zzz,feval(hfun2_80,zzz),'Linewidth',1.6');hold on;plot(zzz,feval(hfun2_100,zzz),'Linewidth',1.6');
% plot(zzz,feval(rfun2_80,zzz),'Linewidth',1.6');plot(zzz,feval(rfun2_100,zzz),'Linewidth',1.6');
% xlabel('z',FS,14);
% leg = legend('$h2, T=80$','$h2, T=100 $',...
% 		'$r2, T=80$','$r2, T=100 $' );
% set(leg,'Interpreter','latex');
% set(leg,FS,14, 'Location','best');
% %D = f1.PaperPosition; % Returns 1x4 vector [left right width height]
% %f1.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
% % print -f1 -dpsc2 sets_drift1.eps;
% print(gcf,'-depsc2','-loose','sets_alpha_z_960_h2.eps');
% %print('-depsc2', 'sets_drift1.eps')


% FS = 'fontsize';
% set(gcf,'paperpositionmode','auto')
% figure;plot(zzz,feval(hfun1_80,zzz).*sigma_y(1)+feval(hfun2_80,zzz).*sigma_y(2)...
% 	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6');hold on;
% plot(zzz,feval(hfun1_100,zzz).*sigma_y(1)+feval(hfun2_100,zzz).*sigma_y(2)...
% 	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6');
% plot(zzz,feval(onefun,zzz).*(mu_y_hat+beta_hat.*feval(zfun,zzz)),'Linewidth',1.6');
% xlabel('z',FS,14);
% lab=ylabel('$\sigma_y \cdot H +\hat{\alpha_y}+\hat{\beta} z$',FS,14);
%  set(lab,'Interpreter','latex');
% legend1 = legend(' T=80','T=100','baseline');
% set(legend1,FS,14, 'Location','best');
% print(gcf,'-depsc2','-loose','sets_alpha_z_960_drift1.eps');


% FS = 'fontsize';
% set(gcf,'paperpositionmode','auto')
% figure;plot(zzz,feval(hfun1_80,zzz).*sigma_z(1)+feval(hfun2_80,zzz).*sigma_z(2)...
% 	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6');hold on;
% plot(zzz,feval(hfun1_100,zzz).*sigma_z(1)+feval(hfun2_100,zzz).*sigma_z(2)...
% 	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6');
% plot(zzz,feval(onefun,zzz).*(mu_z_hat-kappa_hat.*feval(zfun,zzz)),'Linewidth',1.6');
% xlabel('z',FS,14);lab=ylabel('$\sigma_z \cdot H+\hat{\alpha_z}-\hat{\kappa} z$',FS,14);
%  set(lab,'Interpreter','latex');
% legend1 = legend(' T=80','T=100','baseline');
% set(legend1,FS,14, 'Location','best');
% print(gcf,'-depsc2','-loose','sets_alpha_z_960_drift2.eps');


% figure;plot(theta_space(1:kk-1),log(2)./max(rho1(:,1:kk-1)) ,'-*', LW,1.6); hold on;  
% plot(theta_space(1:kk-1),log(2)./max(rho2(:,1:kk-1)),'-*', LW,1.6); 
% xlabel('\theta', FS, 14); ylabel('Half-life',FS, 14);h=legend('H*','R*','Location', 'best'); set(h,FS, 14); 
% print(gcf, '-depsc2', '-loose','sets_hl_theta_alpha_z_960_2.eps')


% figure;plot(theta_space,log(2)./max(rho1), LW,1.6); hold on;  plot(theta_space,log(2)./max(rho2),'r', LW,1.6);
% xlabel('\theta', FS, 14); ylabel('Half-life',FS, 14);h=legend('H*','R*','Location', 'best'); set(h,FS, 14);
% print(gcf, '-depsc2', '-loose','sets_hl_theta_alpha_z_960_2.eps')
% save 'hl_theta_alpha_z_960_3.mat'

%save 'hl_theta_alpha_z_960_2.mat'
% xgrid = linspace(-1,1,5000);
% onefun = chebfun(@(x) 1, dom_z);
% mu = feval(mu_z_hat.*onefun,xgrid)-feval(kappa_hat.*zfun,xgrid);
% sigma = norm(sigma_z).*onefun(xgrid);

% dxc = [xgrid(2)-xgrid(1) 0.5*(xgrid(3:end)-xgrid(1:end-2)) xgrid(end)-xgrid(end-1)];
% dsigma = [sigma(2)-sigma(1) 0.5*(sigma(3:end)-sigma(1:end-2)) sigma(end)-sigma(end-1)]./dxc;

% dlogq = (mu - dsigma.*sigma)*2 ./ (sigma.^2) .* dxc;
% logq = [0 0.5*cumsum(dlogq(1:end-1)+dlogq(2:end))];
% logq = logq - max(logq) + 1; % normalization
% % q = exp(logq)/sum(exp(logq).*dxc - 0.5*exp(logq(1))*dxc(1) - 0.5*exp(logq(end))*dxc(end));
% q = exp(logq)/sum(0.5*sum(exp(logq(1:end-1)).*dxc(1:end-1) + exp(logq(2:end)).*dxc(2:end)));

% figure; 
% plot(xgrid,q, 'Linewidth',1.6);
% xlabel('z',FS,14);
% print('-depsc2','-loose', 'sets_density_alpha_z(0).eps')

% xgrid = linspace(-1,1,5000);
% onefun = chebfun(@(x) 1, dom_z);
% mu = feval(f(z_sol),xgrid)-feval(z_sol.*zfun,xgrid);
% sigma = norm(sigma_z).*onefun(xgrid);

% dxc = [xgrid(2)-xgrid(1) 0.5*(xgrid(3:end)-xgrid(1:end-2)) xgrid(end)-xgrid(end-1)];
% dsigma = [sigma(2)-sigma(1) 0.5*(sigma(3:end)-sigma(1:end-2)) sigma(end)-sigma(end-1)]./dxc;

% dlogq = (mu - dsigma.*sigma)*2 ./ (sigma.^2) .* dxc;
% logq = [0 0.5*cumsum(dlogq(1:end-1)+dlogq(2:end))];
% logq = logq - max(logq) + 1; % normalization
% % q = exp(logq)/sum(exp(logq).*dxc - 0.5*exp(logq(1))*dxc(1) - 0.5*exp(logq(end))*dxc(end));
% q = exp(logq)/sum(0.5*sum(exp(logq(1:end-1)).*dxc(1:end-1) + exp(logq(2:end)).*dxc(2:end)));

% figure; 
% plot(xgrid,q, 'Linewidth',1.6);%hold on; plot(xgrid,q1, 'Linewidth',1.6)
% xlabel('z',FS,14);%legend('distorted','baseline');
% print('-depsc2','-loose', 'sets_density_alpha_z.eps')

% dom_z = [-1 1];
% valfun = chebfun(@(z) z.^3, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
% %valfun= w{3000};
% zfun = chebfun(@(z) z, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
% u = 0.01*feval(diff(f,1),zz)-diff(valfun,1).*zfun;
% z_roots = roots(u);
% roots_fun = polyfit(zz,z_roots,nrank,domain(d));
% kappa_fun = inv(roots_fun);
% z_sol = chebfun({@(z) min(zz),@(z) feval(kappa_fun,z), @(z) max(zz) }, [-1 min(roots_fun) max(roots_fun) 1], ...
% 		'splitting', 'on','eps', 1e-6, 'vectorize')
%plot alpha_y(kappa)
%figure; plot(zz,feval(f,zz))
%plot alpha_y(z)
%figure; plot(z_sol)
%plot d alpha_y(z)/ d kappa
%figure; plot(diff(f,1))