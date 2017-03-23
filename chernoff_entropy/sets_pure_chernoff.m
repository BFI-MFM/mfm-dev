clear all; close all;

pc = parcluster('local');
pc.JobStorageLocation =  strcat('/home/vzhorin/', getenv('SLURM_JOB_ID'));
%pc.JobStorageLocation =  strcat('/home/vzhorin/matlabjobs);
matlabpool(pc, 16)


mu_y_hat = 0.386;
mu_z_hat = 0;
beta_hat = 1;
kappa_hat = 0.019;

mu_y_set = [.25];
mu_set_length = length(mu_y_set);
%mu_y_set = [.25 .27 .30 .33];
%mu_y = mu_y_hat;
for kk = 1:mu_set_length 
	mu_y = mu_y_set(kk);

	beta = beta_hat;
	sigma_y = [0.488 0]';
	sigma_z = [0.013 0.028]';
	sigma = [sigma_y';sigma_z'];

	n_kappa = 500; n_mu = 500;
	half_life_time = 160;

	kappa = linspace(0.0001,0.08,n_kappa);
	mu_z = linspace(-0.012,0,n_mu);
	inv_s = inv(sigma);
	d1 = [0 1];
	xx = [];yy_min = []; yy_max = [];
	for ii = 1:n_kappa
		uu = kappa_hat-kappa(ii);
		icount = 1;yy = [];
		parfor jj=1:n_mu
			mm =  mu_z(jj)-mu_z_hat;
			u = [beta-beta_hat; uu];
			mu = [mu_y-mu_y_hat;mm];
			aa = inv(sigma)*u;
			bb = inv(sigma)*mu;

			zeta2 =  chebfun(@(s) (s-s.^2).*norm(aa).^2,d1,'eps',1e-4,'splitting', 'on','vectorize');
		    zeta1  =  chebfun(@(s) (s-s.^2).*dot(aa,bb), d1,'eps',1e-4,'splitting', 'on','vectorize');
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
	            yy(icount) = mu_z(jj);
	            icount = icount+1;
	        end
	    end
		if ~isempty(yy)
	        xx(ii) = kappa(ii); yy_min(ii) = min(yy); yy_max(ii) = max(yy);
	    end
		disp(ii);

	end

%xxx(kk,:) = xx(:); yyy_min(kk,:) = yy_min(:);
end
FS = 'fontsize';


figure; 
plot(xx, yy_min, 'Linewidth',1.6);hold on; 
%plot(kappa_hat-xxx(2,:), yyy_min(2,:), 'Linewidth',1.6);
%plot(kappa_hat-xxx(3,:), yyy_min(3,:), 'Linewidth',1.6);
xlabel('\kappa',FS,14); ylabel('\mu_z',FS,14);
%legend('\mu_y = 0.33','\mu_y = 0.35','\mu_y = 0.38');
%legend('\mu_y = 0.25','\mu_y = 0.27','\mu_y = 0.30','\mu_y = 0.33');
%set(legend,FS,14, 'Location','best');

print('-depsc2', 'sets_pure_chernoff_mu1.eps');



 figure; 
 plot(xx1, yy1, 'Linewidth',1.6);hold on; 
 plot(xx2, yy2, 'Linewidth',1.6); plot(xx3, yy3, 'Linewidth',1.6);
 plot(xx4, yy4, 'Linewidth',1.6); plot(xx5, yy5, 'Linewidth',1.6);
% %plot(kappa_hat-xxx(2,:), yyy_min(2,:), 'Linewidth',1.6);
% %plot(kappa_hat-xxx(3,:), yyy_min(3,:), 'Linewidth',1.6);
 xlabel('\kappa',FS,14); ylabel('\mu_z',FS,14);
 legend('\mu_y = 0.30','\mu_y = 0.32','\mu_y = 0.34','\mu_y = 0.36','\mu_y = 0.38');
 %legend('\mu_y = 0.25','\mu_y = 0.27','\mu_y = 0.30','\mu_y = 0.33');
 set(legend,FS,14, 'Location','best');

 print('-depsc2', 'sets_pure_chernoff_mu.eps');

save pure_chernoff_160_mu1.mat
matlabpool close