function omega = scale(x,x0,xx,dx,muX,sigmaX,param)

xx_gt_x0 = xx(xx>x0);
omega_integral_x_gt_x0 = cumsum( 2.*muX(xx_gt_x0) ./ sigmaX(xx_gt_x0).^2 .* dx);
omega_x_gt_x0 = @(x) exp( - ...
    interp1(xx_gt_x0, omega_integral_x_gt_x0, x) ...
);

xx_lt_x0 = fliplr( xx(xx<x0) );
omega_integral_x_lt_x0 = cumsum( 2.*muX(xx_lt_x0) ./ sigmaX(xx_lt_x0).^2 .* dx);
omega_x_lt_x0 = @(x) exp(  ...
    interp1(xx_lt_x0, omega_integral_x_lt_x0, x) ...
);

omega2 = omega_x_gt_x0(x(x>x0));
omega1 = omega_x_lt_x0(x(x<x0));

omega = [omega1 omega2];