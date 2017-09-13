function [den,cumden] = Density(process,xx,cutoff,param,imposing0,modelname)
% The file is trying to solve the Fokker-Plance ODE:
%        0 = - d/dx [muX h(x)] + 0.5*d^2/dx^2 [sigmaX^2 h(x)]
% which is equivalent to
%        J* = [muX h(x)] - 0.5*d^2/dx^2 [sigmaX^2 h(x)]
% If imposing0 is empty, the code will let J* = 0
% If imposing0 is not empty, the code will shoot h(imposing0) = 0 by varying J*.

%% Step 1 Read the input

muX = @(x) process.muX(x);
sigmaX = @(x) process.sigmaX(x);

dx = xx(2) - xx(1);

%% Step 2 Solve left boundary problem

xxr = xx(xx>cutoff);
xxl = xx(xx<cutoff);

if (isempty(imposing0))
    Jstar = 0;
else
    if imposing0 > cutoff
        tozero = @(Jstar) FKODE(Jstar,muX,sigmaX,xxr(1),imposing0);
    else
        tozero = @(Jstar) FKODE(Jstar,muX,sigmaX,xxl(end),imposing0);
    end
    Jstar = fsolve(tozero,0.5);
end    
    
toint = @(x,y) muX(x) ./ sigmaX(x).^2 * 2 * y - Jstar;
if imposing0 > cutoff
    tstr = ode23s(toint,[xxr(1) xxr(end)],1);
    tstl = ode23s(toint,[xxl(end) xxl(1)],tstr.y(1)/sigmaX(xxr(1)).^2.*sigmaX(xxl(end)).^2);
else
    tstl = ode23s(toint,[xxl(end) xxl(1)],1);
    tstr = ode23s(toint,[xxr(1) xxr(end)],tstl.y(1)/sigmaX(xxl(end)).^2.*sigmaX(xxr(1)).^2);
end

tmpHr = deval(xxr,tstr);
tmpHl = deval(xxl,tstl);

tmpH = ([tmpHl tmpHr]);
xx = xx(xx ~= cutoff);
tmpdown = sigmaX(xx).^2;
    
% figure;plot(xx,tmpH/max(tmpH)*max(tmpdown));
% hold on;plot(xx,tmpdown);

tmph = 2*tmpH ./ tmpdown;
tmph(tmph<0) = 0;

xx = xx(~isnan(tmph)) ;
tmph = tmph(~isnan(tmph));

den.x = xx;
den.d = tmph ./ sum(tmph*dx);

cumden.x = xx;
cumden.cd = cumsum(den.d*dx);

disp(['Mean of the state is ',num2str( Expectation(den,@(x) x) )]);

