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
%sigma_y = 0.01.*[0.488 0]';
%sigma_z = 0.01.*[0.013 0.028]';
sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';
%alpha = 0.01; 
delta = 0.002; 

sigma = [sigma_y'; sigma_z'];
%ss = [0;0.01];
ss = [0; 0.01];
%ss_set=[0.005 0.01 0.01];
%ss = [0;1];
alpha = norm(inv(sigma)*ss);
ns = 60;
chernoff_s = linspace(0.5,0.9,ns);
%chernoff_s = linspace(0.8,0.95,ns);

n_theta = 60;
theta_space = linspace(10,500,n_theta);

for kk=1:n_theta
	disp('computing theta');
		%kk=1;
	theta = theta_space(kk);
	disp(theta);
	d = [-1 0 1];
	N = chebop(d);
	zfun = chebfun('z', d,'eps', 1e-9, 'vectorize');
	v = chebfun('z.^3', d,'eps', 1e-3, 'vectorize');

	ksi = chebfun({@(z) -alpha*z,@(z) alpha*z}, d, ...
		'splitting', 'on','eps', 1e-6, 'vectorize');

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
		xx = linspace(d(1),d(3),12);
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
    
    file_save = ['sets_density-theta-xi-0-01-' num2str(kk) '.mat'];
	%save(file_save);

	N_c = chebop(d);N_cr = chebop(d);
	%e = chebfun('z.^3', d,'eps', 1e-3, 'vectorize');
	% N.op = 	@(e) -0.5.*chernoff_s.*(1-chernoff_s).*(hfun1_20(v_20).^2+...
	% 	2.*hfun1_20(v_20).*hfun2_20(v_20)+hfun2_20(v_20).^2) +...
	% 	(chernoff_s.*(sigma_z(1).*hfun1_20(v_20)+sigma_z(2).*hfun2_20(v_20))).*diff(e,1)+...
	% 	1/2.*norm(sigma_z).^2.*diff(e,2)+1/2.*norm(sigma_z).^2.*diff(e,1).^2;


	for jj = 1:ns
		s = chernoff_s(jj);
		tic;
		N_c.op = 	@(e) -0.5.*s.*(1-s).*e.*(hfun1(v).^2+...
		+hfun2(v).^2) +...
		(s.*(sigma_z(1).*hfun1(v)+sigma_z(2).*hfun2(v))).*diff(e,1)+...
		1/2.*norm(sigma_z).^2.*diff(e,2);

    	N_c.bc = 'dirichlet';
 		[V,D]= eigs(N_c,100);
 		zz =  -sort(diag(D),'descend');
 		disp(zz(1));
 		rho(jj,kk) = real(zz(1));

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
end 	
	save(file_save);
	rho_max = max(rho);
	theta_fun = chebfun.spline(theta_space(1:30),rho_max);
	save ('sets-scan_theta-xi-0-01.mat', 'theta_fun', 'rho_max', 'theta_space')
 	FS = 'fontsize';
 	set(gcf,'paperpositionmode','auto')
 	figure;
 	 plot(theta_space(1:30), log(2)./rho_max, 'Linewidth',1.6);
	xlabel('\theta',FS,14);ylabel('T(\theta)',FS,14);
	print(gcf, '-depsc2', '-loose','sets_chernoff_findtheta-xi-01.eps')


	% figure;
 % 	plot(log(2)./theta_fun1, 'Linewidth',1.6);hold on;
 % 	plot(log(2)./theta_fun2, 'Linewidth',1.6);
 % 	plot(log(2)./theta_fun3, 'Linewidth',1.6);
 % 	plot(linspace(0,200,100), 80.*ones(1,100),'--');
	% xlabel('\theta',FS,14);ylabel('T(\theta)',FS,14);
	 legend1 = legend('\xi=0.005', '\xi=0.01');
	 set(legend1,FS,14, 'Location','best');

	 print(gcf, '-depsc2', '-loose','sets_chernoff_findtheta-xi.eps')






