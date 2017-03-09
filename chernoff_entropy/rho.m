function [rho] = rho(s, param, param_hat)
%rho= rhos(s) Summary of this function goes here
%   Detailed explanation goes here

mu_y_hat = param_hat(1);
mu_z_hat = param_hat(2);
beta_hat = param_hat(3);
kappa_hat = param_hat(4);

mu_y = param(1);
mu_z = param(2);
beta = param(3);
kappa = param(4);


sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';
sigma = [sigma_y';sigma_z'];
inv_s = inv(sigma);
uu = kappa_hat-kappa;
u = [beta-beta_hat; uu];
mm = mu_z-mu_z_hat;
mu = [mu_y-mu_y_hat;mm];
aa = inv_s*u;
bb =  inv_s*mu;
        
zeta2 = (s-s.^2).*norm(aa).^2;
zeta1 =  (s-s.^2).*dot(aa,bb);
zeta0 =  (s-s.^2).*norm(bb).^2;
%disp('zeta computed');disp(num2str(zeta0));
dd = sqrt((kappa_hat-s.*uu).^2+zeta2.*norm(sigma_z).^2);
lambda2 = (kappa_hat-s.*uu-dd)./norm(sigma_z).^2;
%disp('lambda2 computed');
lambda1 =  (mu_z_hat*(1-s)+s.*mu_z.*lambda2-zeta1)./dd;
%figure; plot(lambda1);
%disp('lambda1 computed');
rho =  zeta0./2 -1/2.*norm(sigma_z).^2.*lambda2-...
    1/2.*norm(sigma_z).^2.*lambda1.^2-(mu_z_hat*(1-s)+s.*mu_z).*lambda1;

end

