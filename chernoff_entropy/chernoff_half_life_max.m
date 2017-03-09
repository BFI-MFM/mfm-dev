function delta = chernoff_half_life_max(param, param_hat, hl)
%rho= rhos(s) Summary of this function goes here
%   Detailed explanation goes here

% mu_y_hat = 0.386;
% mu_z_hat = 0;
% beta_hat = 1;
% kappa_hat = 0.019;

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

n_s = 10000;
ss = linspace(0,1,n_s);
inv_s = inv(sigma);
xx = []; yy_min = []; yy_max = [];
ii=1;
uu(ii) = kappa_hat-kappa;
jj=1;
mm(jj)= mu_z-mu_z_hat;
u =[beta-beta_hat; uu(ii)];
    icount = 1;yy = [];
        %mu= [0; mm(jj)];
        mu = [mu_y-mu_y_hat;mm(jj)];
        %aa = inv_s(2,2)*u*0.01;
        aa = inv_s*u;
        bb =  inv_s*mu;
        rho = zeros(n_s,1);

        for kk=1:n_s
            s =ss(kk);
            zeta2 = (s-s.^2).*norm(aa).^2;
            zeta1 =  (s-s.^2).*dot(aa,bb);
            zeta0 =  (s-s.^2).*norm(bb).^2;
            %disp('zeta computed');disp(num2str(zeta0));
            dd = sqrt((kappa_hat-s.*uu(ii)).^2+zeta2.*norm(sigma_z).^2);
            lambda2 = (kappa_hat-s.*uu(ii)-dd)./norm(sigma_z).^2;
            %disp('lambda2 computed');
            lambda1 =  (mu_z_hat*(1-s)+s.*mu_z(jj).*lambda2-zeta1)./dd;
            %figure; plot(lambda1);
            %disp('lambda1 computed');
            rho(kk) =  zeta0./2 -1/2.*norm(sigma_z).^2.*lambda2-...
                1/2.*norm(sigma_z).^2.*lambda1.^2-(mu_z_hat*(1-s)+s.*mu_z(jj)).*lambda1;
            %disp('rho computed');
        end
        
         rr = max(rho);
        %disp(log(2.)/rr);
        t = log(2)./rr;
        delta = abs(t - hl);