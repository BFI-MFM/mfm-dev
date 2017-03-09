clear all;close all;
set_path_tensor
FS = 'fontsize';
load hl_theta_beta_960_hl-60-1

xgrid = linspace(-1,1,5000);
onefun = chebfun(@(x) 1, dom_z);
mu = mu_z-kappa.*zfun(xgrid);
sigma = norm(sigma_z).*onefun(xgrid);

dxc = [xgrid(2)-xgrid(1) 0.5*(xgrid(3:end)-xgrid(1:end-2)) xgrid(end)-xgrid(end-1)];
dsigma = [sigma(2)-sigma(1) 0.5*(sigma(3:end)-sigma(1:end-2)) sigma(end)-sigma(end-1)]./dxc;

dlogq = (mu - dsigma.*sigma)*2 ./ (sigma.^2) .* dxc;
logq = [0 0.5*cumsum(dlogq(1:end-1)+dlogq(2:end))];
logq = logq - max(logq) + 1; % normalization
% q = exp(logq)/sum(exp(logq).*dxc - 0.5*exp(logq(1))*dxc(1) - 0.5*exp(logq(end))*dxc(end));
q = exp(logq)/sum(0.5*sum(exp(logq(1:end-1)).*dxc(1:end-1) + exp(logq(2:end)).*dxc(2:end)));

figure; 
plot(xgrid,q, 'Linewidth',1.6);
xlabel('z',FS,14);
%print('-depsc2','-loose', 'sets_beta_hl_290_z_density.eps');

qfun =  chebfun.spline(xgrid,q);	

d1 = [ hfun1(dom_z(1)) hfun1(-0.00001)];
d2 = [ hfun2(dom_z(1)) hfun2(-0.00001)];

xx1 = linspace(d1(1),d1(end),10000);
xx2 = linspace(d2(1),d2(end),10000);

zh11 = inv(chebfun(@(x) feval(hfun1,x), [dom_z(1) -0.00001], 'eps', 1e-6, 'splitting','on'),'eps', 1e-3,'splitting','on');
disp('zh1 computed');
zh21 = inv(chebfun(@(x) feval(hfun2,x),  [dom_z(1) -0.00001], 'eps', 1e-6, 'splitting','on'),'eps', 1e-3,'splitting','on');
disp('zh2 computed');
g_infty11 = feval(qfun,zh11(xx1)).*abs(feval(diff(zh11,1),xx1));
g_infty21 = feval(qfun,zh21(xx2)).*abs(feval(diff(zh21,1),xx2));

figure; %plot(xx1,g_infty1,'.','Linewidth',1.6);hold on; 
h2=area(xx1,g_infty11);
h2.FaceColor = [0 0.5 0.5];
h2.EdgeColor = 'none';
%plot(xx2,g_infty2,'.','Linewidth',1.6);
legend1 = legend('p(h_1)');
set(legend1,FS,14, 'Location','best');
xlabel('h1',FS,14);
print('-depsc2', 'sets_beta_hl_60_h_density.eps')

