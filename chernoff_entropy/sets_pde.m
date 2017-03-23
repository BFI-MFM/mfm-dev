clear all; close all;
tic;
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
ss = [0;0.01];
%ss = [0;1];
alpha = norm(inv(sigma)*ss);
 % pc = parcluster('local');
 % matlabpool(pc, 3)

% mu = [0.386 0];beta = 1;
% k = 0.019; sigma = [0.013 0.028];
% alpha = 0.01; delta = 0.002;
%cheboppref.setDefaults('display','off')
%theta = 100;
%d = [-inf inf];
theta_space = [1000000];

for kk=1:1
	theta = theta_space(kk);
	d = [-1 0 1];
	N = chebop(d);
	zfun = chebfun('z', d,'eps', 1e-9, 'vectorize');
	%norm_f = @(v) alpha*zfun.*sqrt((sigma(1)+sigma(2).*diff(v,1)).^2) 
	%v = chebfun('(z-2).^3', d,'eps', 1e-3, 'vectorize');
	v = chebfun('z.^3', d,'eps', 1e-3, 'vectorize');
	%v = chebfun({@(z) -z.^2,@(z) z.^2}, d,'splitting', 'on','eps', 1e-9, 'vectorize')

	ksi = chebfun({@(z) -alpha*z,@(z) alpha*z}, d, ...
		'splitting', 'on','eps', 1e-6, 'vectorize');

	ksi_term = @(v) ksi.*sqrt((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2);
	theta_term = @(v) ((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2)*1./(2.*theta);
	%N.init = chebfun(@(z) exp(-z.^2), d);
	% N.op = @(v) -delta.*v + (mu(1)+beta.*zfun)+diff(v,1).*(mu(2)-k.*zfun)+...
	% 	1/2.*abs(sigma(2)).^2.*diff(v,2)-1./(2.*theta).*...
	% 	(sigma(1).^2+2.*sigma(1).*sigma(2).*diff(v,1)+sigma(2)^.2*diff(v,1).^2)-...
	% 	ksi.*(sigma(1)+sigma(2).*diff(v,1));
	%N1.op =	@(v) -delta.*v + (mu(1)+zfun)-diff(v,1).*k.*zfun+...
	%	1/2.*sigma(2).^2.*diff(v,2)-1./(2.*theta).*(sigma(1).^2+2.*sigma(1).*sigma(2).*diff(v,1));
	N.op = 	@(v) -delta.*v + (mu_y+zfun)-diff(v,1).*k.*zfun+...
		1/2.*norm(sigma_z).^2.*diff(v,2)-ksi_term(v) - theta_term(v);
		%-1./(2.*theta).*sigma(2)^.2*diff(v,1).^2;%-ksi.*(sigma(1)+sigma(2).*diff(v,1));
	w{1} = v;
	for i = 2:50000
		ff = N(w{i-1});
		v =   w{i-1} + ff;
		xx = linspace(-1,1,12);
		%w{i} = v;
		w{i} = chebfun.spline(xx,v(xx));
		disp(i); disp(abs(max(ff))+abs(min(ff)));
	end
	value_fun{kk} = v;
end

hfun0 = @(v) -(1./theta+ksi./sqrt((sigma(1,1)+sigma(2,1)*diff(v,1)).^2+...
		(sigma(2,1)+sigma(2,2)*diff(v,1)).^2));
hfun1 = @(v)	hfun0(v).*(sigma(1,1)+sigma(2,1)*diff(v,1));
hfun2 = @(v)	hfun0(v).*(sigma(2,1)+sigma(2,2)*diff(v,1));

figure;plot(hfun1(v),'r', 'Linewidth',1.6);
hold on;plot(hfun2(v),'b', 'Linewidth',1.6)
% Create xlabel
xlabel('z','FontSize',16);
legend('h_1(z)','h_2(z)','Location','Best');
print('-depsc2', 'sets_h-theta-5.eps')

hhfun1 =  @(v)sigma(1,1).*hfun1(v)+sigma(1,2).*hfun2(v)
hhfun2 =  @(v)sigma(2,1).*hfun1(v)+sigma(2,2).*hfun2(v)
% plot(hhfun1(v),'r','Linewidth',1.6);hold on;plot(hhfun2(v),'b','Linewidth',1.6)
% xlabel('z','FontSize',16);
% legend('\sigma_y(1)*h_1+\sigma_z(1)*h_2','\sigma_y(1)*h_1+\sigma_z(2)*h_2');

 figure;plot(diff(v,1),'Linewidth',1.6)
 xlabel('z','FontSize',16);
 ylabel('v^\prime(z)','FontSize',16);
print('-depsc2', 'sets_upsilon-theta-5.eps')

onefun = chebfun(@(x) 1, d);
figure;plot(sigma_y(1)+sigma_z(1).*diff(v,1),'b','Linewidth',1.6);hold on;
plot(sigma_y(2)+sigma_z(2).*diff(v,1),'r','Linewidth',1.6);
xlabel('z','FontSize',16);
legend('\sigma_y(1)+\sigma_z(1)*v^\prime(z)','\sigma_y(2)+\sigma_z(2)*v^\prime(z)','Location','Best');
print('-depsc2', 'sets_sigmas-theta-5.eps')

figure;plot(mu_y+zfun,'r','Linewidth',1.6');hold on;
plot(mu_y+zfun+hhfun1(v),'b','Linewidth',1.6)
xlabel('z','FontSize',16);
legend('\mu_y+z','\mu_y+z+ \sigma_y(1)*h_1+\sigma_z(1)*h_2','Location','Best');
print('-depsc2', 'sets_diag1-theta-5.eps')

figure;plot(mu_z-k.*zfun,'r','Linewidth',1.6);hold on;
plot(mu_z-k.*zfun+hhfun2(v),'b','Linewidth',1.6)
xlabel('z','FontSize',16);
legend('\mu_z-\kappa z', '\mu_z-\kappa z+\sigma_y(2)*h_1+\sigma_z(2)*h_2','Location','Best');
print('-depsc2', 'sets_diag2-theta-5.eps')
save 'sets-theta-5.mat'
%KFE
% xgrid = linspace(-1,1,10000);
% onefun = chebfun(@(x) 1, d);
% mu = -k.*zfun(xgrid);
% sigma = norm(sigma_z).*onefun(xgrid);

% dxc = [xgrid(2)-xgrid(1) 0.5*(xgrid(3:end)-xgrid(1:end-2)) xgrid(end)-xgrid(end-1)];
% dsigma = [sigma(2)-sigma(1) 0.5*(sigma(3:end)-sigma(1:end-2)) sigma(end)-sigma(end-1)]./dxc;

% dlogq = (mu - dsigma.*sigma)*2 ./ (sigma.^2) .* dxc;
% logq = [0 0.5*cumsum(dlogq(1:end-1)+dlogq(2:end))];
% logq = logq - max(logq) + 1; % normalization
% % q = exp(logq)/sum(exp(logq).*dxc - 0.5*exp(logq(1))*dxc(1) - 0.5*exp(logq(end))*dxc(end));
% q = exp(logq)/sum(0.5*sum(exp(logq(1:end-1)).*dxc(1:end-1) + exp(logq(2:end)).*dxc(2:end)));

% figure; plot(xgrid, q,'Linewidth',1.6) 
% xlabel('Z','FontSize',16);
% print('-depsc2', 'sets_density.eps')