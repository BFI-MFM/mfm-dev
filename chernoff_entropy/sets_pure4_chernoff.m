clear all; close all;
%set_path_tensor
%pc=parcluster('local');
%matlabpool(pc, 12);
mu_y = 0.386;
%mu_y = 0.35;
mu_z =  0;
%mu_z =  -0.002;
beta = 1;
k = 0.019;
%k = 0.01;

sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';

sigma = [sigma_y';sigma_z'];
n_kappa = 400; n_mu = 400;n_s =5000;
uu = k-linspace(0.0,0.1,n_kappa);
mm = linspace(-0.02,0.02,n_mu);
ss = linspace(0,1,n_s);
inv_s = inv(sigma);
%bb = norm(inv_s*mu2);

rho=zeros(n_s,1);
xx = []; yy_min = []; yy_max = [];
for ii=1:n_kappa
    u =[0; uu(ii)];
    icount = 1;yy = [];
    for jj = 1:n_mu
        mu= [0; mm(jj)];
        %mu = [0;-0.002];
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
            dd = sqrt((k-s.*uu(ii)).^2+zeta2.*norm(sigma_z).^2);
            lambda2 = (k-s.*uu(ii)-dd)./norm(sigma_z).^2;
            %disp('lambda2 computed');
            lambda1 =  (s.*mm(jj).*lambda2-zeta1)./dd;
            %figure; plot(lambda1);
            %disp('lambda1 computed');
            rho(kk) =  zeta0./2 -1/2.*norm(sigma_z).^2.*lambda2-1/2.*norm(sigma_z).^2.*lambda1.^2-(s.*mm(jj)).*lambda1;
            %disp('rho computed');
        end
        rr = max(rho);
        %disp('log(2.)/rr');
        if(log(2.)/rr >= 160);
            yy(icount) = mm(jj);
            icount = icount+1;
        end
    end
    if ~isempty(yy)
        xx(ii) = uu(ii); yy_min(ii) = min(yy); yy_max(ii) = max(yy);
    end
    disp(ii);
end

figure;
plot(k-xx, yy_min,'-r');hold on;  plot(k-xx, yy_max,'-b'); xlabel('\kappa'); ylabel('\mu');
print('-depsc2', 'sets_pure_chernoff.eps');
save pure_chernoff_medium.mat


