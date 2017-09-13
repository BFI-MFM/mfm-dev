function BoundaryType(process,domainx,x0,param)

%% Step 1 Read the input

muX = process.muX;
sigmaX = process.sigmaX;

xmin = domainx(1);
xmax = domainx(2);
% xx = linspace(domainx(1),domainx(2),nx);
% dx = xx(2) - xx(1);
% xx = xx(xx~=x0);

%% Step 2 Define scale function

% scale_fcn = @(x) OMEGA(x,x0,xx,dx,muX,sigmaX,param);
pull = @(x) 2.*muX(x) ./ sigmaX(x).^2;
scale_fcn = @(x) exp(-integral(pull,x0,x));

%% Step 3 Integrals near boundaries
    
% I_xl_inner_itgl = cumsum( scale_fcn(xx) .* dx );
% I_xl_inner_all = @(y) interp1(xx,I_xl_inner_itgl,y) .* 2./ sigmaX(y).^2 ./ scale_fcn(y);
% % I_xmintmp = cumsum( I_xl_inner_all( xx(xx<x0) ) .* dx );
% % I_xmintmp = I_xmintmp(~isnan(I_xmintmp));
% % I_xmin = I_xmintmp(end);
% I_xmin = sum( I_xl_inner_all( xx(xx<x0) ) .* dx );
% 
% Jinner = @(x) JInner(x,x0,xx,dx,muX,sigmaX,param,scale_fcn);
% J_inner_all = @(y) Jinner(y) .* 2./ sigmaX(y).^2 ./ scale_fcn(y);
% % J_xmintmp = cumsum( J_inner_all( xx(xx<x0) ) .* dx );
% % J_xmintmp = J_xmintmp(~isnan(J_xmintmp));
% % J_xmin = J_xmintmp(end);
% J_xmin = sum( J_inner_all( xx(xx<x0) ) .* dx );
% 
% I_xu_inner_itgl = cumsum( scale_fcn(fliplr(xx)) .* dx );
% I_xu_inner_all = @(y)  interp1(fliplr(xx),I_xu_inner_itgl,y) .* 2./ sigmaX(y).^2 ./ scale_fcn(y);
% % I_xmaxtmp = cumsum( I_xu_inner_all( xx(xx>x0) ) .* dx );
% % I_xmaxtmp = I_xmaxtmp(~isnan(I_xmaxtmp));
% % I_xmax = I_xmaxtmp(end);
% I_xmax = cumsum( I_xu_inner_all( xx(xx>x0) ) .* dx );
% 
% % J_xmaxtmp = - cumsum( J_inner_all( xx(xx>x0) ) .* dx );
% % J_xmaxtmp = J_xmaxtmp(~isnan(J_xmaxtmp));
% % J_xmax = J_xmaxtmp(end);
% J_xmax = - sum( J_inner_all( xx(xx>x0) ) .* dx );

I_xmin_inner = @(y) integral(scale_fcn,xmin,y,'ArrayValued',true);
I_tmp = @(y) 2*I_xmin_inner(y)./(sigmaX(y).^2.*scale_fcn(y));
I_xmin = integral(I_tmp,xmin,x0,'ArrayValued',true);
I_xmax_inner = @(y) integral(scale_fcn,y,xmax,'ArrayValued',true);
I_tmp = @(y) 2*I_xmax_inner(y)./(sigmaX(y).^2.*scale_fcn(y));
I_xmax = integral(I_tmp,x0,xmax,'ArrayValued',true);

J_xmin_inner = @(y) integral(scale_fcn,y,x0,'ArrayValued',true);
J_tmp = @(y) 2*J_xmin_inner(y)./(sigmaX(y).^2.*scale_fcn(y));
J_xmin = integral(J_tmp,xmin,x0,'ArrayValued',true);
J_xmax_inner = @(y) integral(scale_fcn,y,x0,'ArrayValued',true);
J_tmp = @(y) 2*J_xmax_inner(y)./(sigmaX(y).^2.*scale_fcn(y));
J_xmax = integral(J_tmp,xmax,x0,'ArrayValued',true);

%% Step 4 Classify boundaries

proxy_inf = 999999;
if isnan(I_xmin)
    I_xmin = Inf;
end
if isnan(I_xmax)
    I_xmax = Inf;
end
if isnan(J_xmin)
    J_xmin = Inf;
end
if isnan(J_xmax)
    J_xmax = Inf;
end

if (I_xmin <= proxy_inf && J_xmin <= proxy_inf)
    disp('left boundary is regular');
elseif (I_xmin > proxy_inf && J_xmin <= proxy_inf)
    disp('left boundary is entrance');
elseif (I_xmin <= proxy_inf && J_xmin > proxy_inf)
    disp('left boundary is exit');
elseif (I_xmin > proxy_inf && J_xmin > proxy_inf)
    disp('left boundary is natural');
end

if (I_xmax <= proxy_inf && J_xmax <= proxy_inf)
    disp('right boundary is regular');
elseif (I_xmax > proxy_inf && J_xmax <= proxy_inf)
    disp('right boundary is entrance');
elseif (I_xmax <= proxy_inf && J_xmax > proxy_inf)
    disp('right boundary is exit');
elseif (I_xmax > proxy_inf && J_xmax > proxy_inf)
    disp('right boundary is natural');
end
