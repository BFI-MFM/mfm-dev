format('compact');
clear;
global thetass rhoh m f sigmasqr gama g alpha lambda rho l

g=0.02;rhoh=0.04;rho=rhoh;%rho and rhoh are the same
l=1.84;gama=2; sigma=0.09;m=4;lambda=0.6;

sigmasqr=sigma^2;
f=0.0; % f is proportional management fee paid by households. The final version we set f=0, but we investigated a version with positive f.   
sigmahatsqr=sigmasqr;

aaa=rho+g*(gama-1)+0.5*sigmasqr*gama*(1-gama)-gama*rhoh*l/(1+l);
thetass=m/(1-lambda+m);
aa=rhoh+g*(gama-1)+0.5*sigmasqr*gama*(1-gama);
bb=f+rhoh-rho-g*(gama-1)-0.5*sigmasqr*gama*(1-gama);

T=10^4; yending=(1+l)/rhoh;
F=zeros(T,1);y=zeros(T,1);Fprime=zeros(T,1); F2prime=zeros(T,1);ggg=zeros(T,1);G=zeros(T,1);thetas=zeros(T,1);
riskpremium=zeros(T,1); miuy=zeros(T,1); liquid=zeros(T,1);sigmav=zeros(T,1);corre=zeros(T,1);sharpratio=zeros(T,1);
portiondele=zeros(T,1);psd=zeros(T,1);eta=zeros(T,1);cIvol=zeros(T,1);vol=zeros(T,1);leverage=zeros(T,1);
thetab=zeros(T,1); r=zeros(T,1);r2=zeros(T,1);sigmay=zeros(T,1);eta2=zeros(T,1);res=zeros(T,1);res2=zeros(T,1);alphai=zeros(T,1);

xx=zeros(T,1);xprime=zeros(T,1);x2prime=zeros(T,1);miux=zeros(T,1);sigmax=zeros(T,1);
dy=yending/T;
y=dy*[1:1:T]';
start=floor(0.25/dy); % one can adjust starting value 0.25 for better precision

%% start shooting for the right slope Fprimenew at y0 so that F(yb) lands at (1+l)/rhoh

Fprimett=0;Fprimettt=Fprimett-2;

y0=start*dy;Fprimenew=(Fprimett+Fprimettt)/2;
F0=(1+Fprimenew*l)/(rho+g*(gama-1)+0.5*sigmasqr*gama*(1-gama)-gama*rhoh*l/(1+l));
tol=0.001;
F0try=F0+y0*Fprimenew;
le=0.001;
yb=(1+l)/rhoh-le;Fb=yb; 

yspan=[y0:dy:yb];
x=size(yspan);
ylength=x(2);
options1=odeset('reltol',1e-9);

[Y,H]=ode15s(@delefully2hw15s,yspan,[F0try;Fprimenew],options1);
x=size(Y);
length=x(1);
mis=(ylength==length)*(H(length,1)-Fb)+(ylength>length)*(-100)

while abs(mis)>tol;
if mis<0
    Fprimettt=Fprimenew;
    Fprimenew=0.5*Fprimettt+0.5*Fprimett;
    else 
    Fprimett=Fprimenew;
    Fprimenew=0.5*Fprimettt+0.5*Fprimett;
end;
F0=(1+Fprimenew*l)/(rho+g*(gama-1)+0.5*sigmasqr*gama*(1-gama)-gama*rhoh*l/(1+l));
F0try=F0+y0*Fprimenew;
[Y,H]=ode15s(@delefully2hw15s,yspan,[F0try;Fprimenew],options1);
x=size(Y);length=x(1);mis=(ylength==length)*(H(length,1)-Fb)+(ylength>length)*(-100)
end

%% given the solution, calculate cutoff point for capital constraint yc, risk premium, interest rate, etc 

F(1:start)=F0+Fprimenew.*[0:dy:(start-1)*dy];
Fprime(1:start)=Fprimenew.*ones(start,1);
for t=1:length
    F(start+t-1)=H(t,1);Fprime(start+t-1)=H(t,2);
end
% compute yc 
tt=T;ttt=1;
newtt=floor((tt+ttt)/2);
err=F(newtt)-y(newtt)/thetass;
while tt-ttt>1
    if err>0;
        ttt=newtt;newtt=floor((newtt+tt)/2);
    else
        tt=newtt;newtt=floor((newtt+ttt)/2);
    end
    err=F(newtt)-y(newtt)/thetass;    
end
tc=tt;delta=(F(tt)-F(ttt))/(y(tt)-y(ttt));
yc=(F(tt)-y(ttt)*delta)/(1/thetass-delta);
Fc=yc/thetass;
% calculate risk premium and interest rate etc
totallength=length-1+start;
for x=totallength+1:T
    F(x)=F(totallength)+(x-totallength)*dy;Fprime(x)=1;
end
for t=2:T-1
    thetas(t)=(y(t)<=thetass*F(t))*(1-lambda)*y(t)/(F(t)-lambda*y(t))+(y(t)>thetass*F(t))*m/(1+m);
    F2prime(t)=(Fprime(t)-Fprime(t-1))/dy;
    thetab(t)=y(t)-thetas(t)*F(t);
    ggg(t)=1/(1-thetas(t)*Fprime(t));Gprime=thetas(t)*ggg(t)^2*(Fprime(t)-Fprime(t-1))/dy;
    riskpremium(t)=gama*sigmasqr*(1+rhoh*thetab(t)*ggg(t)/(1+l-rhoh*y(t)))...
        *(1-Fprime(t)*ggg(t)*thetab(t)/F(t));
    r(t)=((-m*f*(y(t)>thetass*F(t))+rho+g*gama)*(1+l-rhoh*y(t))-rhoh*gama*ggg(t)*(thetas(t)+l-thetab(t)*g-rhoh*y(t)-f*min((1-lambda)*y(t),m*(F(t)-y(t)))+0.5*sigmasqr*thetab(t)^2*Gprime)...
    -0.5*gama*(gama+1)*sigmasqr*(1+l-rhoh*y(t))*(1+rhoh*thetab(t)*ggg(t)/(1+l-rhoh*y(t)))^2)/(1+l-rhoh*y(t)+rhoh*gama*ggg(t)*thetab(t));
    sigmay(t)=-thetab(t)*sigma*ggg(t);
    miuy(t)=ggg(t)*(thetas(t)+l+(r(t)+sigmasqr-g)*thetab(t)-rhoh*y(t)+0.5*sigmasqr*Gprime*thetab(t)^2-f*min((1-lambda)*y(t),m*(F(t)-y(t))));
    
    %transform y to x
    xx(t)=1-y(t)/F(t);
    xprime(t)=(-F(t)+y(t)*Fprime(t))/(F(t))^2;
    x2prime(t)=(y(t)*F2prime(t)*(F(t))^2+2*F(t)^2*Fprime(t)-2*y(t)*F(t)*(Fprime(t))^2)/(F(t)^4);
    miux(t)=xprime(t)*miuy(t)+(1/2)*x2prime(t)*sigmay(t)^2;
    sigmax(t)=xprime(t)*sigmay(t);
    
    r2(t)=-m*f*(y(t)>thetass*F(t))+rho+gama*(g-rhoh*(miuy(t)+sigmay(t)*sigma)/(1+l-rhoh*y(t)))-gama*(gama+1)/2*(sigma-rhoh*sigmay(t)/(1+l-rhoh*y(t)))^2;
    
    cIvol(t)=sigma-rhoh*sigmay(t)/(1+l-rhoh*y(t));
    portiondele(t)=min((1-lambda)*y(t),m*(F(t)-y(t)))/F(t);
    liquid(t)=gama*(cIvol(t)^2-sigmasqr);
    sigmav(t)=sigma*(1-Fprime(t)*thetab(t)/F(t)*ggg(t));
    sharpratio(t)=riskpremium(t)/sigmav(t);
    corre(t)=1/sqrt(1+sigmahatsqr/sigmav(t)^2);
    psd(t)=(1+l-rhoh*yc)^gama/(1+l-rhoh*y(t))^gama;
    eta(t)=(sigma*(F(t)-y(t))+(Fprime(t)-1)*sigmay(t))*F(t)/(F(t)-y(t))/(sigma*F(t)+Fprime(t)*sigmay(t));
    eta2(t)=1/(1-(m+1)*thetab(t)/F(t));
    res(t)=(1+l-rho*F(t)+(Fprime(t)-1)*miuy(t)+0.5*F2prime(t)*sigmay(t)^2+...
         (F(t)-y(t))*((1-gama)*g+0.5*gama*(gama-1)*sigmasqr)+(Fprime(t)-1)*sigma*sigmay(t)*(1-gama))*(1+l-rhoh*y(t))+...
         gama*rhoh*miuy(t)*(F(t)-y(t))+0.5*gama*(1+gama)*rhoh^2*sigmay(t)^2*(F(t)-y(t))/(1+l-rhoh*y(t))+...
         (Fprime(t)-1)*gama*sigmay(t)^2*rhoh+gama*rhoh*sigmay(t)*(1-gama)*sigma*(F(t)-y(t));
    M=gama*rhoh*miuy(t)/(1+l-rho*y(t))+0.5*gama*(1+gama)*rhoh^2*sigmay(t)^2/(1+l-rhoh*y(t))^2+(1-gama)*g+0.5*gama*(gama-1)*sigmasqr+...
         gama*rhoh*sigmay(t)*(1-gama)*sigma/(1+l-rho*y(t));
    N=miuy(t)+gama*sigmay(t)^2*rhoh/(1+l-rho*y(t))+sigma*sigmay(t)*(1-gama);
    res2(t)=l-y(t)*M-N;
    alphai(t)=1+thetab(t)/(F(t)-y(t));
    leverage(t)=1-1/alphai(t);
end
alphai(T)=alphai(T-1);alphai(1)=1;
x=100;
figure(1)
subplot(421)
plot(xx(2:totallength-x),F(2:totallength-x))
xlabel('w/P: scaled specialist wealth x')
title('equilibrium Price/Dividend ratio F')
subplot(422)
plot(xx(2:totallength-x),r2(2:totallength-x))
title('interest rate')
xlabel('w/P: scaled specialist wealth x')
subplot(423)
plot(xx(2:totallength-x),riskpremium(2:totallength-x))
title('risk premium')
xlabel('w/P: scaled specialist wealth x')
subplot(424)
plot(xx(2:totallength-x),sigmav(2:totallength-x))
title('stock volatility')
xlabel('w/P: scaled specialist wealth x')
subplot(425)
plot(xx(2:totallength-x),sharpratio(2:totallength-x))
title('Sharpe ratio')
xlabel('w/P: scaled specialist wealth x')
subplot(426)
plot(xx(2:totallength-x),alphai(2:totallength-x));
xlabel('w/P: scaled specialist wealth x')
title('Intermediary�s Position in Risky Asset')
subplot(427)
plot(xx(2:totallength-500),miux(2:totallength-500));
xlabel('w/P: scaled specialist wealth x')
title('Drift of x: \mu_x')
subplot(428)
plot(xx(2:totallength-500),sigmax(2:totallength-500));
xlabel('w/P: scaled specialist wealth x')
title('Diffusion of x: \sigma_x')

save('baselinesolution.mat');


