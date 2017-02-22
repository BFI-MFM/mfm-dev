%solving F'' using G; used in implicitschemefully2.m 
function dh=delefully2hw15s(y,h)
global thetass rhoh m f sigmasqr gama g  lambda rho l
dh=zeros(2,1);
%y=1.5;h=[2;1];hprime=[2;1];
%g=0.001;rhoh=0.02; lambda=0.1;
%f=0.03; gama=2; sigmasqr=0.01;mhat=4;
if y<=thetass*h(1)
    thetas=(1-lambda)*y/(h(1)-lambda*y);
else
  % thetass=m*(1-y/h(1));
  thetas=m/(1+m);
 
end

ggg=1/(1-thetas*h(2));

thetab=y-thetas*h(1);
xxx=(1+l-rhoh*y)*(ggg-1)/thetas/h(1)+gama*rhoh*ggg;
yyy=thetas+l+thetab*(g*(gama-1)+rho-m*f*(y>thetass*h(1)))-rhoh*y-f*min((1-lambda)*y,m*(h(1)-y));
zzz=1+rhoh*ggg*thetab/(1+l-rhoh*y);
fff=0.5*thetab^2*sigmasqr*ggg*(1+l+rhoh*y*(gama-1))/h(1)/(1+l-rhoh*y+rhoh*gama*ggg*thetab)*ggg^2;
hhh=rho+g*(gama-1)-1/h(1)-m*f*(y>thetass*h(1));
iii=0.5*gama*(1-gama)*sigmasqr*zzz*(y-thetab*ggg)*(1+l-rhoh*y-rhoh*ggg*thetab)/(1+l-rhoh*y+rhoh*gama*ggg*thetab)/thetas/h(1);    
jjj=-xxx*yyy/(1+l-rhoh*y+rhoh*gama*ggg*thetab);
%kkk=(1+l-rho*y)*(thetab*ggg-y)*m*f*(y>thetass*h(1))/(thetas*h(1)*(1+rho*gama*ggg*thetab));
% 
% if y>(1+l)/rhoh-0.02;
%     zzz
% end
% 
%  
dh(1)=h(2);
dh(2)=(hhh+iii+jjj)/fff; % don't put fff in front of hprime(2) 
%ff(2)=hprime(2)*fff-(hhh+iii+jjj);



  