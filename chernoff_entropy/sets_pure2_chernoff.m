clear all; close all;

% create a local cluster object
pc = parcluster('local');

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
matlabpool(pc,16)


mu_y_hat = 0.386;
mu_z_hat = 0;
beta_hat = 1;
kappa_hat = 0.019;

mu_y = mu_y_hat;
%mu_y=0.35;

beta = beta_hat;
sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';
sigma = [sigma_y';sigma_z'];

n_kappa = 500; n_mu = 5000;
half_life_time = 120;

kappa = linspace(0.0001,0.07,n_kappa);
mu_z = linspace(-0.008,0.008,n_mu);
inv_s = inv(sigma);
d1 = [0 1];
xx = [];yy_min = []; yy_max = [];
for ii = 1:n_kappa
	icount = 1;yy = [];
	uu = kappa_hat-kappa(ii);
	parfor jj=1:n_mu
		mm =  mu_z(jj)-mu_z_hat;
		u = [beta-beta_hat; uu];
		mu = [mu_y-mu_y_hat;mm];
		aa = inv(sigma)*u;
		bb = inv(sigma)*mu;

		zeta2 =  chebfun(@(s) (s-s.^2).*norm(aa).^2,d1,'eps',1e-4,'splitting', 'on','vectorize');
	    zeta1 =  chebfun(@(s) (s-s.^2).*dot(aa,bb), d1,'eps',1e-4,'splitting', 'on','vectorize');
		zeta0 =  chebfun(@(s) (s-s.^2).*norm(bb).^2 ,d1,'eps',1e-4,'splitting', 'on','vectorize');
		dd =  chebfun(@(s) sqrt((kappa_hat-s.*uu).^2+feval(zeta2,s).*norm(sigma_z).^2),...
			d1,'eps',1e-4,'splitting', 'on','vectorize');
		%disp('zeta computed');
		lambda2 =  chebfun(@(s) (kappa_hat-s.*uu-feval(dd,s))./norm(sigma_z).^2,...
			d1, 'eps',1e-4,'splitting', 'on','vectorize');
		%disp('lambda2 computed');
		lambda1 =  chebfun(@(s)(mu_z_hat.*(1-s)+...
			s.*mu_z(jj).*feval(lambda2,s)-feval(zeta1,s))./feval(dd,s),...
			d1,'eps',1e-4,'splitting', 'on','vectorize');
		%figure; plot(lambda1);
		%disp('lambda1 computed');
		rho =  chebfun(@(s)  feval(zeta0,s)./2 -1/2.*norm(sigma_z).^2.*feval(lambda2,s)-...
				1/2.*norm(sigma_z).^2.*feval(lambda1,s).^2- (mu_z_hat.*(1-s)+s.*mu_z(jj)).*feval(lambda1,s),...
			d1,'eps',1e-4,'splitting', 'on','vectorize');
		%disp('rho computed');
		rr(ii,jj) = max(rho);
		disp(jj);
	% if(log(2.)/rr(ii,jj) >= half_life_time);
  	%           yy(icount) = mm;
  	%           icount = icount+1;
  	%       end
	end
	for jj=1:n_mu
		if(log(2.)/rr(ii,jj) >= half_life_time);
            yy(icount) =mu_z(jj);
            icount = icount+1;
        end
    end
	if ~isempty(yy)
        xx(ii) = uu; yy_min(ii) = min(yy); yy_max(ii) = max(yy);
    end
	disp(ii);

end

FS = 'fontsize';
load('victor.mat')

figure; plot(bound(:,1),bound(:,3),'-k'); hold on;plot(bound(:,1),bound(:,2),'-m');
plot(kappa_hat-xx, yy_max,'-b', 'Linewidth',1.6);hold on;  plot(kappa_hat-xx, yy_min,'-r','Linewidth',1.6);
xlabel('\kappa',FS,14); ylabel('\mu',FS,14);
legend('upper bound(Lloyd)','lower bound(Lloyd)','upper bound(Victor)','lower bound(Victor)')
set(legend,FS,14, 'Location','best');

print('-depsc2', 'sets_pure_chernoff_160_2.eps');

save pure_chernoff_160_2.mat
matlabpool close