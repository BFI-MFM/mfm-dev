clear all; close all;
tic;
set_path_tensor
%load sets_density_new.mat
% setup parameters
%mu_y = .0386;
mu_y = 3.86;
mu_z =  0;
beta = 1;
k = 0.019; 

sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';

delta = 0.002; 

sigma = [sigma_y'; sigma_z'];
%ss = [0; 0.01];
%zz=linspace(-100,100,100);
%ss =[-0.52+0.038.*zz; -0.0055+0.006.*zz]
%alpha = norm(inv(sigma)*ss);

ns = 50;
chernoff_s = linspace(0.5,0.999,ns);
n_theta = 15;
theta_space = linspace(25,80,n_theta);

for kk=1:n_theta
	disp('computing theta');
		%kk=1;
	theta = theta_space(kk);
	disp(theta);
	d = [-3 3];
	N = chebop(d);
	zfun = chebfun('z', d,'splitting', 'on','eps', 1e-6, 'vectorize');
	v = chebfun('z', d,'splitting', 'on','eps', 1e-6, 'vectorize');
	ss =[-0.052+0.038.*zfun; -0.0055+0.006.*zfun];
	ff = inv(sigma)*ss;
	
	ksi = sqrt(ff{1}.^2+ff{2}.^2);
	%ksi = chebfun({@(z) -alpha*z,@(z) alpha*z}, d, ...
	%	'splitting', 'on','eps', 1e-6, 'vectorize');

	ksi_term = @(v) ksi.*sqrt((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2);
	theta_term = @(v) ((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2)*1./(2.*theta);
	N.op = 	@(v) -delta.*v + (mu_y+zfun)-diff(v,1).*k.*zfun+...
		1/2.*norm(sigma_z).^2.*diff(v,2)-ksi_term(v) - theta_term(v);
	w{1} = v;
	for i = 2:3000
		ff = N(w{i-1});
		v =   w{i-1} + ff;
		xx = linspace(d(1),d(2),24);
		%w{i} = v;
		w{i} = chebfun.spline(xx,v(xx));
		%disp(i); disp((abs(max(ff))+abs(min(ff)))./abs(max(v)));
	end
	disp(i); disp((abs(max(ff))+abs(min(ff)))./abs(max(v)));
	hfun0 = @(v) -(1./theta+ksi./sqrt((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2));
	hfun1 = @(v)	hfun0(v).*(sigma(1,1)+sigma(2,1)*diff(v,1));
	hfun2 = @(v)	hfun0(v).*(sigma(2,1)+sigma(2,2)*diff(v,1));
	
	rfun0 = @(v) -(ksi./sqrt((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2));

	rfun1 = @(v)	rfun0(v).*(sigma(1,1)+sigma(2,1)*diff(v,1));
	rfun2 = @(v)	rfun0(v).*(sigma(2,1)+sigma(2,2)*diff(v,1));
    
    file_save = ['sets_density-theta-xi-second-part2-' num2str(kk) '.mat'];


	N_c = chebop(d);N_cr = chebop(d);

	for jj = 1:ns
		s = chernoff_s(jj);
		tic;
		N_c.op = 	@(e) -0.5.*s.*(1-s).*e.*(hfun1(v).^2+...
		+hfun2(v).^2) +...
		(s.*(sigma_z(1).*hfun1(v)+sigma_z(2).*hfun2(v))).*diff(e,1)+...
		1/2.*norm(sigma_z).^2.*diff(e,2);

    	N_c.bc = 'neumann';
 		[V,D]= eigs(N_c,100);
 		zz =  -sort(real(diag(D)),'descend');
 		disp(zz(1));
 		rho(jj,kk) = zz(1);

 	% 	N_cr.op = 	@(e) -0.5.*s.*(1-s).*e.*(rfun1(v).^2+...
		% +rfun2(v).^2) +...
		% (s.*(sigma_z(1).*rfun1(v)+sigma_z(2).*rfun2(v))).*diff(e,1)+...
		% 1/2.*norm(sigma_z).^2.*diff(e,2);

  %   	N_cr.bc = 'neumann';
 	% 	[V,D]= eigs(N_cr,100);
 	% 	zz = -sort(diag(D),'descend');
 	% 	lambda(jj,kk) = zz(1);


% 		ee{jj} = V{1};
 		disp(jj);
 		toc;
 	end
 	save(file_save);
end 	
	save(file_save);
	rho_max = max(rho);
	theta_fun = chebfun.spline(theta_space,rho_max);
	save ('sets-scan_theta-xi-second-part2.mat', 'theta_fun', 'rho_max', 'theta_space')
 	FS = 'fontsize';
%  	set(gcf,'paperpositionmode','auto')
%  	figure;
%  	 plot(theta_space, log(2)./rho_max, 'Linewidth',1.6);
% 	xlabel('\theta',FS,14);ylabel('T(\theta)',FS,14);
% 	print(gcf, '-depsc2', '-loose','sets_chernoff_findtheta-xi-second.eps')

% FS = 'fontsize';
%  	set(gcf,'paperpositionmode','auto')
%  	figure;
%  	 plot(chernoff_s, rho1, 'Linewidth',1.6);
% xlabel('s',FS,14);
% legend1 = legend('\theta=35', '\theta=\infty');
% set(legend1,FS,14, 'Location','best');
% 	print(gcf, '-depsc2', '-loose','sets_chernoff_rho-xi_theta-all-second-80.eps')


	% figure;
 % 	plot(log(2)./theta_fun1, 'Linewidth',1.6);hold on;
 % 	plot(log(2)./theta_fun2, 'Linewidth',1.6);
 % 	plot(log(2)./theta_fun3, 'Linewidth',1.6);
 % 	plot(linspace(0,200,100), 80.*ones(1,100),'--');
	% xlabel('\theta',FS,14);ylabel('T(\theta)',FS,14);
	% legend1 = legend('\xi=0.005', '\xi=0.01', '\xi=0.015');
	% set(legend1,FS,14, 'Location','best');

	% print(gcf, '-depsc2', '-loose','sets_chernoff_findtheta-xi.eps')






