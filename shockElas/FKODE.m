function out = FKODE(Jstar,muX,sigmaX,cutoff,imposing0)

toint = @(x,y) muX(x) ./ sigmaX(x).^2 * 2 .* y - Jstar;

tst = ode23s(toint,[cutoff imposing0],1);
out = tst.y(end) ./ sigmaX(tst.x(end)).^2;


