clear all;close all;
mu_y_hat = 0.386;
mu_z_hat = 0;
beta_hat = 1;
kappa_hat = 0.019;
param_hat = [mu_y_hat  mu_z_hat beta_hat kappa_hat];

n_mu_y = 1000; mu_y_low = 0.001;
mu_y_domain =linspace(mu_y_low, mu_y_hat, n_mu_y);

n_mu_z = 1; mu_z_low =mu_z_hat;mu_z_high= mu_z_hat;
mu_z_domain = linspace(mu_z_low, mu_z_high, n_mu_z);

n_beta = 1; beta_high =beta_hat; beta_low = beta_hat;
beta_domain =linspace(beta_low, beta_high, n_beta);

n_kappa = 1; kappa_low = kappa_hat;
kappa_domain =linspace(kappa_low, kappa_hat, n_kappa);
half_life = 120;
ii=[];jj=[];kk=[];ll=[];
for i = 1:n_mu_y
    %mu_y = 0.3;
    param(1) = mu_y_domain(i);
    tic;
    for j = 1:n_mu_z
        %mu_z = 0;
        param(2) = mu_z_domain(j);
        for k = 1:n_beta
            % beta = 1;
            param(3) = beta_domain(k);
            for l = 1:n_kappa
                %kappa = 0.019;
                param(4) = kappa_domain(l);
                t(l) = chernoff_half_life_max(param, param_hat);
            end
            ww = [];icount = 0;
            for l = 1:n_kappa
                if t(l) >= half_life
                    icount =icount+1;
                    ww(icount) = kappa_domain(l);
                end
            end
            if ~isempty(ww)
                xx(i) =mu_y_domain(i); ww_min(i) = min(ww);
                param(1) =   xx(i); param(4) =   ww_min(i);
                tt(i) = chernoff_half_life_max(param, param_hat);
                disp(xx(i)); disp(ww_min(i)); disp(tt(i));
            end
        end
    end
    toc;
    disp(i);
end

