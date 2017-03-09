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
ss = [0; 0.005];
%ss_set=[0.005 0.01 0.015];
ss_set=[0.005 0.01];
%ss = [0;1];
alpha = norm(inv(sigma)*ss);
ns = 100;
chernoff_s = linspace(0.001,0.999,ns);
%chernoff_s = linspace(0.8,0.95,ns);

%theta_space = [18.43 38 109.14];
theta_space = [14.3 10000000 27.88 10000000];
for ii=1:2
	ss(2) = ss_set(ii);
	alpha = norm(inv(sigma)*ss);
	theta = theta_space(ii);
for kk=1:2
	%kk=1;
	theta = theta_space((2*ii-1)+(kk-1));
	disp((2*ii-1)+(kk-1));
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
	disp(ii); disp((abs(max(ff))+abs(min(ff)))./abs(max(v)));
	hfun0 = @(v) -(1./theta+ksi./sqrt((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2));
	hfun1 = @(v)	hfun0(v).*(sigma(1,1)+sigma(2,1)*diff(v,1));
	hfun2 = @(v)	hfun0(v).*(sigma(2,1)+sigma(2,2)*diff(v,1));
	
	rfun0 = @(v) -(ksi./sqrt((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2));

	rfun1 = @(v)	rfun0(v).*(sigma(1,1)+sigma(2,1)*diff(v,1));
	rfun2 = @(v)	rfun0(v).*(sigma(2,1)+sigma(2,2)*diff(v,1));
    
	N_c = chebop(d);N_cr = chebop(d);
	%e = chebfun('z.^3', d,'eps', 1e-3, 'vectorize');
	% N.op = 	@(e) -0.5.*chernoff_s.*(1-chernoff_s).*(hfun1_20(v_20).^2+...
	% 	2.*hfun1_20(v_20).*hfun2_20(v_20)+hfun2_20(v_20).^2) +...
	% 	(chernoff_s.*(sigma_z(1).*hfun1_20(v_20)+sigma_z(2).*hfun2_20(v_20))).*diff(e,1)+...
	% 	1/2.*norm(sigma_z).^2.*diff(e,2)+1/2.*norm(sigma_z).^2.*diff(e,1).^2;


	for jj = 1:ns
		s = chernoff_s(jj);
%		tic;
		N_c.op = 	@(e) -0.5.*s.*(1-s).*e.*(hfun1(v).^2+...
		+hfun2(v).^2) +...
		(s.*(sigma_z(1).*hfun1(v)+sigma_z(2).*hfun2(v))).*diff(e,1)+...
		1/2.*norm(sigma_z).^2.*diff(e,2);

    	N_c.bc = 'neumann';
 		[V,D]= eigs(N_c,100);
 		zz =  -sort(diag(D),'descend');
 		disp(zz(1));
 		rho(jj,ii) = zz(1);

 	% 	N_cr.op = 	@(e) -0.5.*s.*(1-s).*e.*(rfun1(v).^2+...
		% +rfun2(v).^2) +...
		% (s.*(sigma_z(1).*rfun1(v)+sigma_z(2).*rfun2(v))).*diff(e,1)+...
		% 1/2.*norm(sigma_z).^2.*diff(e,2);

  %   	N_cr.bc = 'neumann';
 	% 	[V,D]= eigs(N_cr,100);
 	% 	zz = -sort(diag(D),'descend');
 	% 	lambda(jj,ii) = zz(1);


% 		ee{jj} = V{1};
 		%disp(jj);
% 		toc;
 	end
 	file_save = ['sets_density-theta-xi-all2-' num2str(ii) '-'  num2str(kk) '.mat'];
	save(file_save);

end 	
end

 	FS = 'fontsize';
 	set(gcf,'paperpositionmode','auto')
 	figure;
 	plot(chernoff_s, rho1(:,1),'Linewidth',1.6);hold on;plot(chernoff_s, rho1(:,2),'Linewidth',1.6);
 	plot(chernoff_s, rho1(:,3),'Linewidth',1.6);plot(chernoff_s, rho1(:,4),'Linewidth',1.6);
	xlabel('s',FS,14);ylabel('\rho(s)',FS,14);
	legend1 = legend('\xi=0.005, \theta=14.3', ...
		'\xi=0.005, \theta=\infty','\xi=0.01,\theta=27.9', '\xi=0.01,\theta=\infty');
	set(legend1,FS,14, 'Location','best');
	print(gcf, '-depsc2', '-loose','sets_chernoff_rho-xi_theta-all2.eps')


	% figure;
 % 	plot(chernoff_s, lambda(:,1),'Linewidth',1.6);hold on;plot(chernoff_s, lambda(:,2),'Linewidth',1.6);
 % 	plot(chernoff_s, lambda(:,3),'Linewidth',1.6);
	% xlabel('s',FS,14);ylabel('\lambda(s)',FS,14);
	% legend1 = legend('\xi=0.005', '\xi=0.01', '\xi=0.015');
	% set(legend2,FS,14, 'Location','best');
	% print(gcf, '-depsc2', '-loose','sets_chernoff_lambda-xi.eps')




