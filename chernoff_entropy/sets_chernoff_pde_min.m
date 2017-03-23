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

% load kappa_chernoff_set_120.mat;
% ff = squeeze(ff)';
% zz(ff==max(ff))=[];
% ff(ff==max(ff))=[];
% ff1 =ff;
% zz1 =zz;
% load kappa_chernoff_set_240.mat;
% ff = squeeze(ff)';
% zz(ff==max(ff))=[];
% ff(ff==max(ff))=[];
% ff2 =ff;
% zz2 =zz;
% load kappa_chernoff_set_360.mat;
% ff = squeeze(ff)';
% zz(ff==max(ff))=[];
% ff(ff==max(ff))=[];
% ff3 =ff;
% zz3 =zz;
%  figure; plot(zz1,ff1,LW,1.6); hold on; 
%  plot(zz2,ff2,LW,1.6); plot(zz3,ff3,LW,1.6);
% xlabel('\kappa', FS, 14); ylabel('\alpha_y',FS, 14); 
%  legend('T=120', 'T=240', 'T=360');

load kappa_chernoff_set_360.mat;
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
nrank = 24;
f= polyfit(zz, ff,nrank,domain(d));
plot(zz, feval(f,zz)-ff,'r',LW,1.6);

n_theta = 1;
%theta_space = linspace(80,100,n_theta);
theta_space = 100000;
ss_sigma  = sigma*sigma';
dom_z = [-1 0 1];
onefun = chebfun(@(z) 1, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');

for kk=1:n_theta
	disp('computing theta');
		%kk=1;
	theta = theta_space(kk);
	disp(theta);
	valfun = chebfun(@(z) z+cos(0.3.*z-1), dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
	%valfun= w{3000};
	zfun = chebfun(@(z) z, dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
	w{1} = valfun;
	for i = 2:200
		u = 0.01.*feval(diff(f,1),zz)-diff(valfun,1).*zfun;
		z_roots = roots(u);
		roots_fun = polyfit(zz,z_roots,nrank,domain(d));
		kappa_fun = inv(roots_fun);
		z_sol = chebfun({@(z) min(zz),@(z) feval(kappa_fun,z), @(z) max(zz) }, [dom_z(1) min(roots_fun) max(roots_fun) dom_z(end)], ...
			'eps', 1e-6, 'vectorize');
        f.domain = z_sol.domain;
		N = chebop(dom_z);
		% N.op = 	@(v) -delta.*v + 0.01.*(f(z_sol)+beta.*zfun)+diff(v,1).*(mu_z-z_sol.*zfun)+...
		% 	1/2.*norm(sigma_z).^2.*diff(v,2) - theta_term(v);

 		 theta_term = @(v)  (0.01.*(0.01*ss_sigma(1,1)+diff(valfun,1).*ss_sigma(1,2))+...
			 diff(valfun,1).*(0.01.*ss_sigma(2,1)+diff(valfun,1).*ss_sigma(2,2)))...
			 *1./(2.*theta);
		N.op = 	@(v) -delta.*v + 0.01.*(f(z_sol)+beta.*zfun)+diff(v,1).*(mu_z-z_sol.*zfun)+...
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

	[feval(diff(v,1), -1) 0.01/(delta+min(z_sol))]
	[feval(diff(v,1), 1) 0.01/(delta+max(z_sol))]
	%[feval(diff(v,1), 20) 0.01/(delta+max(z_sol))]

	figure; plot(f(z_sol).*0.01-z_sol.*zfun, LW,1.6); 
	xlabel('z', FS, 14); ylabel('0.01*\alpha_y-z\kappa(z)',FS, 14); 
	print(gcf, '-depsc2', '-loose','sets_drift1.eps')

	figure; plot(f(z_sol)+zfun, LW,1.6); 
	xlabel('z', FS, 14); ylabel('\alpha_y+z',FS, 14); 
	print(gcf, '-depsc2', '-loose','sets_drift2.eps')

		figure; plot(diff(valfun,1), LW,1.6); 
	xlabel('z', FS, 14); ylabel('d v/ dz',FS, 14); 
	print(gcf, '-depsc2', '-loose','sets_drift_val_derivative.eps')

	figure; plot(-z_sol.*zfun, LW,1.6);hold on;
	 plot(-kappa_hat.*zfun, LW,1.6);
	  xlabel('z', FS, 14); 
	leg=legend('$-\kappa(z) z$','$-\hat{\kappa} z$', 'Location', 'best'); 
	set(leg,FS, 14); set(leg,'Interpreter','latex');
	print(gcf, '-depsc2', '-loose','sets_drift_hat.eps')

		figure; plot(f(z_sol), LW,1.6);hold on;
	 plot(mu_y_hat.*onefun, LW,1.6);
	  xlabel('z', FS, 14); leg=legend('$\alpha_y(z)$','$\hat{\alpha_y}$', 'Location', 'best'); 
	  set(leg,FS, 14); set(leg,'Interpreter','latex');
	print(gcf, '-depsc2', '-loose','sets_drift_hat2.eps')

	rfun1 = sigma_inv(1,1).*(f(z_sol)-mu_y_hat)+sigma_inv(1,2).*(mu_z-mu_z_hat)+...
		sigma_inv(1,1).*(beta-beta_hat).*zfun+sigma_inv(1,2).*(kappa_hat-z_sol).*zfun;
	
	rfun2 = sigma_inv(2,1).*(f(z_sol)-mu_y_hat)+sigma_inv(2,2).*(mu_z-mu_z_hat)+...
		sigma_inv(2,1).*(beta-beta_hat).*zfun+sigma_inv(2,2).*(kappa_hat-z_sol).*zfun;
   
    figure;plot(-rfun1,LW,1.6);hold on;plot(-rfun2,'r',LW,1.6); xlabel('z', FS, 14); ylabel('-r^*(z)',FS, 14); 
    h=legend('-r_1(z)','-r_2(z)','Location', 'best'); set(h,FS, 14); 
    print(gcf, '-depsc2', '-loose','sets_drift_distortion.eps')
    %file_save = ['sets_density-theta-xi-0-01-' num2str(kk) '.mat'];
	%save(file_save);


	parfor jj = 1:ns
		N_c = chebop(dom_z);
		s = chernoff_s(jj);
		tic;
		N_c.op = 	@(e) -0.5.*s.*(1-s).*e.*(rfun1(v).^2+rfun2(v).^2) +...
		(s.*(sigma_z(1).*rfun1(v)+sigma_z(2).*rfun2(v)) + mu_z_hat-kappa_hat.*zfun).*diff(e,1)+...
		1/2.*norm(sigma_z).^2.*diff(e,2);

    	N_c.bc = 'dirichlet';
 		[V,D]= eigs(N_c,100);
 		e_v =  -sort(real(diag(D)),'descend');
 		disp(e_v(1));
 		rho(jj,kk) = e_v(1);
 		disp(jj);
 		toc;
 	end
end 

delete(par_obj)
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