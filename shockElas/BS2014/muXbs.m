function out = muXbs(x,param)

out = x .* ((sigmaXbs(x,param)./x).^2 + (param.a+1/param.phi)./param.Q(x) - 1/param.phi - param.rho);

end