function dq = F(x,q,param)

a = param.a ;
abar = param.abar ;
delta = param.delta ;
phi = param.phi ;
sigma = param.sigma ;
rho = param.rho ;
rhobar = param.rhobar ;

dq = q./(Psi(x,q,param)-x) .* ( 1 - sigma .* ( q.* (Psi(x,q,param)-x)./(x.*(1-x).*(a-abar)) ).^0.5 );

