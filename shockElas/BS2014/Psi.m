function psi = Psi(x,q,param)

a = param.a ;
abar = param.abar ;
delta = param.delta ;
phi = param.phi ;
sigma = param.sigma ;
rho = param.rho ;
rhobar = param.rhobar ;

psi = ( q.*((1-x).*rhobar+x.*rho) + 1/phi.*(q-1) - abar ) ./ (a-abar) ;
end
