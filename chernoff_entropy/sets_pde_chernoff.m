clear all; close all;
tic;

mu_y_hat = .386;
mu_z_hat =  0;
beta_hat = 1;
kappa_hat = 0.019; 

mu_y = mu_y_hat;
mu_y = 0.38;
mu_z = mu_z_hat;
mu_z = 0.05;
beta = beta_hat;
beta = 1.01;
kappa = 0.01;

sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';
delta = 0.002; 

rho_tilde = log(2)/80;
sigma = [sigma_y'; sigma_z'];
inv_s = inv(sigma);
%alpha = norm(inv(sigma)*ss);

mu1 = mu_y-mu_y_hat;
mu2 = mu_z-mu_z_hat;
mu = [mu1;mu2];
u1 = beta-beta_hat;
u2 = kappa_hat-kappa;
u = [u1;u2];
aa = inv_s*mu;
bb =  inv_s*u;
d1 = -10; d2 = 10;
z = linspace(d1,d2,100);
eta = @(zz) aa*ones(1,size(zz,2))+bb*zz;
norm_eta = @(zz) norm(aa)^2+2*dot(aa,bb)*zz+norm(bb)^2*zz.^2;
figure;plot(z,sqrt(norm_eta(z)));
figure;plot(z,eta(z));
%kink detection


            % zeta2 = (s-s.^2).*norm(aa).^2;
            % zeta1 =  (s-s.^2).*dot(aa,bb);
            % zeta0 =  (s-s.^2).*norm(bb).^2;