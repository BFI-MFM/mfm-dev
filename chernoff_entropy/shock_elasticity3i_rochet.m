clear all; close all;
tic;
set_path_tensor
%for elasticity of E
load /home/vzhorin/mfm/rochet/rochet_logE.mat;
%for elasticity of R(E)
%load /home/vzhorin/mfm/rochet/rochet.mat;
global aa_diff bb_diff cc_diff psi delx

% emax = 0.105;
% ee= emax/(d(2)-d(1));
% beta = (xfun-param.p).*Lfun;%./(ee.*(d(2)+0.001-xfun));
% alpha = param.sigma0.*Lfun;%./(ee.*(d(2)+0.001-xfun));

%beta = mucheb.*xi+g -sigma^2/2+1/2*sigmacheb.^2.*diff(xi,1);
%alpha = sigmacheb.*xi+sigma;

%start discretization
dt = 0.05; tmax= 1;
nt = tmax/dt+1;
%tt = logspace(-3,2,nt);
tt = linspace(0,tmax,nt);

low_bound = logspace(-6,-2,5);
grid_size= 1000;%*logspace(4);
el_final =zeros(10,grid_size-1);
yp0 = linspace(0,0,grid_size)';
for kk = 1:5
x0 = domain_E(1)+low_bound(kk);xmax = domain_E(2);
%delx = (xmax-x0)/(nt-1); 
nx = grid_size;
%delx = (log(xmax)-log(x0))/nx;
bfun = mucheb+alpha.*sigmacheb;
sfun = 2*bfun./sigmacheb.^2;
sfun_inv = inv(sfun);
%xx = exp(linspace(log(x0),log(xmax),nx))'; 
xx = linspace(sfun(x0),sfun(xmax),nx)'; 

delx = [xx(2:end)-xx(1:end-1); xx(end)-xx(end-1)];
aa_diff = 1/2.*sigmacheb(xx).^2./delx.^2;
bb_diff = (mucheb(xx)+sigmacheb(xx).*alpha(xx))./delx;
cc_diff = beta(xx)+alpha(xx).^2/2;
psi(1) = xi(x0);
psi(2) = xi(xmax);

u0 = linspace(1,1,nx)'; 
r_r=[0.0234 0.0272 0.0324 ]; 
r= sort(feval(inv(Rfun),r_r),2, 'ascend');
num_r = size(r,2);
for i=1:num_r    
    ir(i)=min(find(xx>=r(i)));
end
for i=1:nt-1	
	tspan = [tt(i) tt(i+1)];
%	yp0 = linspace(0,0,nx)';
	[u0,yp0] = decic('ela_rhs_boundary_fixed',tspan(1),u0',[],yp0,[]);
	[t uu] = ode15i('ela_rhs_boundary_fixed',tspan,u0',yp0);	
	nn = size(t);
	u(i,:) = uu(nn(1),:);
	sh_el(i,1:nx-1) = diff(log(u(i,:)))./diff(xx)'.*sigmacheb(xx(1:nx-1)')+alpha(xx(1:nx-1)');
	u0(:) = u(i,:);
    disp(['time:' num2str(tt(i))]); 
    disp(['state:' num2str(xx(ir)') ' sh el:' num2str(sh_el(i,ir))]);
end 
 disp(['kk:' num2str(kk)]); 
sh_el(nt,:) = sh_el(nt-1,:);
el_final(kk,:) = sh_el(nt,:);
end

for kk = 1:1
	x0 = domain_E(1)+low_bound(kk);xmax = domain_E(2);
	x1(kk,:) = linspace(x0,xmax,nx-1)'; 
end

plot(x1(1,:),el_final(1,:));hold on;
plot(x1(2,:),el_final(2,:));
plot(x1(3,:),el_final(3,:));
plot(x1(4,:),el_final(4,:));
plot(x1(5,:),el_final(5,:));
legend('min(E) = 0.000001','min(E) = 0.00001','min(E) = 0.0001','min(E) = 0.001'...
,'min(E) = 0.01');

semilogy(x1(1,:),el_final(1,:));hold on;
semilogy(x1(2,:),el_final(2,:));
semilogy(x1(3,:),el_final(3,:));
semilogy(x1(4,:),el_final(4,:));
semilogy(x1(5,:),el_final(5,:));
legend('min(E) = 0.000001','min(E) = 0.00001','min(E) = 0.0001','min(E) = 0.001'...
,'min(E) = 0.01');



% surf(tt(1:i-1),xx(1:1900),sh_el(1:i-1,1:1900)'), shading interp, lighting phong, axis tight
%  view([-90 90]), colormap(autumn); set(gca,'zlim',[0.1 0.5])
%  light('color',[1 1 0],'position',[0,0,1])
%  material([0.30 0.60 0.60 40.00 1.00])

% plot(xx(1:nx-1),sh_el(1,1:nx-1),xx(1:nx-1),ones(1,nx-1)*0.1,'--r'), grid on
%  MS = 'markersize';
%  pp=interp1(xx(1:nx-1),sh_el(1,1:nx-1), r,'pcchip');
%  hold on, plot(r,pp,'.r',MS,30);


% plot(xx(1:4999),sh_el(1,1:4999),xx(1:4999),ones(1,4999)*0.1,'--r'), grid on
% MS = 'markersize';
% r=[0.05 0.091 0.27]; pp=interp1(xx(1:4999),sh_el(1,1:4999), r,'pcchip');
% hold on, plot(r,pp,'.r',MS,30);

% plot(xx(1:4999),pr_el(1,1:4999),xx(1:4999),ones(1,4999)*0.2,'--r'), grid on
% MS = 'markersize';
% r=[0.05 0.091 0.27]; pp=interp1(xx(1:4999),pr_el(1,1:4999), r,'pcchip');
% hold on, plot(r,pp,'.r',MS,30);

 efinal = size(sh_el,1)-1;
 for ii=1:num_r
    for jj=1:efinal
    ee(ii,jj) = interp1(xx(1:end-1),sh_el(jj,1:end), r(ii),'pcchip');
%   ep1(jj) = interp1(xx(1:end-1),pr_el(jj,1:end), r(1),'pcchip');
%   ep2(jj) = interp1(xx(1:end-1),pr_el(jj,1:end), r(2),'pcchip');
%   ep3(jj) = interp1(xx(1:end-1),pr_el(jj,1:end), r(3),'pcchip');
end
 end

 plot(tt(1:efinal),ee(1:3,:), 'Linewidth',1.6), grid on;
 legend('25%', '50%','75%')
 print('-depsc2', ['sh_el_rochet.eps'])
 
% MS = 'markersize';
% r=[0.05 0.091 0.27]; pp=interp1(xx(1:9900),pr_el(1,1:9900), r,'pcchip');
% hold on, plot(r,pp,'.r',MS,30);

save rochet-implicit2.mat  