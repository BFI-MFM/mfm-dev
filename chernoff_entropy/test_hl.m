tic;
mu_y_hat = 0.386;
mu_z_hat = 0;
beta_hat = 1;
kappa_hat = 0.019;
param_hat = [mu_y_hat  mu_z_hat beta_hat kappa_hat];

mu_y = 0.3;
mu_z = 0;
beta = 1;
kappa = 0.019;
param = [mu_y  mu_z beta kappa];

 t = chernoff_half_life(param, param_hat)
 toc
