function Rfun = rochet_fun(domain_E, param)
% Solving ODE for the equilibrium loan rate as a function of aggregate equity 

N = chebop(domain_E);
alpha = (param.R-param.p).^(-param.beta);
Lfun = @(r) alpha*(param.R-r).^param.beta;
diff_Lfun = @(r) -alpha.*param.beta.*(param.R-r).^(param.beta-1);

N.op = @(e,r) diff(r) + ...
	1/param.sigma0^2.*(2*(param.rho-param.r)*param.sigma0^2+...
	(r-param.p-param.r).^2+2*param.r*e*(r-param.p-param.r)./Lfun(r))/...
	(Lfun(r)-diff_Lfun(r).*(r-param.p-param.r));
N.rbc = param.p+param.r;
Rfun =  N\domain_E(1);