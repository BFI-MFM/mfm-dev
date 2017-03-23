clear all;close all;
load '/home/vzhorin/new_ss/sets_density-theta-xi-all2-2-1.mat'

v_20 = v;
theta_20 = theta;
sigma_20 = sigma;
ksi_20 =ksi;

hfun0_20 = @(v) -(1./theta_20+ksi_20./sqrt((sigma_20(1,1)+sigma_20(2,1)*diff(v,1)).^2+...
		(sigma_20(2,1)+sigma_20(2,2)*diff(v,1)).^2));
rfun0_20 = @(v) -(ksi_20./sqrt((sigma_20(1,1)+sigma_20(2,1)*diff(v,1)).^2+...
		(sigma_20(2,1)+sigma_20(2,2)*diff(v,1)).^2));

hfun1_20 = @(v)	hfun0_20(v).*(sigma_20(1,1)+sigma_20(2,1)*diff(v,1));
hfun2_20 = @(v)	hfun0_20(v).*(sigma_20(2,1)+sigma_20(2,2)*diff(v,1));

rfun1_20 = @(v)	rfun0_20(v).*(sigma_20(1,1)+sigma_20(2,1)*diff(v,1));
rfun2_20 = @(v)	rfun0_20(v).*(sigma_20(2,1)+sigma_20(2,2)*diff(v,1));

hhfun1_20 =  @(v)sigma_20(1,1).*hfun1_20(v)+sigma_20(1,2).*hfun2_20(v);
hhfun2_20 =  @(v)sigma_20(2,1).*hfun1_20(v)+sigma_20(2,2).*hfun2_20(v);

load '/home/vzhorin/new_ss/sets_density-theta-xi-all2-2-2.mat'
v_infty = v;
theta_infty = theta;
ksi_infty =ksi;
sigma_infty = sigma_20;
hfun0_infty = @(v) -(1./theta_infty+ksi_infty./sqrt((sigma_infty(1,1)+sigma_infty(2,1)*diff(v,1)).^2+...
		(sigma_infty(2,1)+sigma_infty(2,2)*diff(v,1)).^2));
rfun0_infty = @(v) -(ksi_infty./sqrt((sigma_infty(1,1)+sigma_infty(2,1)*diff(v,1)).^2+...
		(sigma_infty(2,1)+sigma_infty(2,2)*diff(v,1)).^2));


hfun1_infty = @(v)	hfun0_infty(v).*(sigma_infty(1,1)+sigma_infty(2,1)*diff(v,1));
hfun2_infty = @(v)	hfun0_infty(v).*(sigma_infty(2,1)+sigma_infty(2,2)*diff(v,1));

rfun1_infty = @(v)	rfun0_infty(v).*(sigma_infty(1,1)+sigma_infty(2,1)*diff(v,1));
rfun2_infty = @(v)	rfun0_infty(v).*(sigma_infty(2,1)+sigma_infty(2,2)*diff(v,1));

hhfun1_infty =  @(v)sigma_infty(1,1).*hfun1_infty(v)+sigma_infty(1,2).*hfun2_infty(v);
hhfun2_infty =  @(v)sigma_infty(2,1).*hfun1_infty(v)+sigma_infty(2,2).*hfun2_infty(v);

FS = 'fontsize';
f1= figure;
set(gcf,'paperpositionmode','auto')
plot(mu_y+zfun,'Linewidth',1.6');hold on;
plot(mu_y+zfun+hhfun1_infty(v_infty),'Linewidth',1.6);
plot(mu_y+zfun+hhfun1_20(v_20),'Linewidth',1.6);
xlabel('z',FS,14);
legend1 = legend('original model','\theta=\infty, \xi=0.01','\theta=27.77, \xi=0.01');
set(legend1,FS,14, 'Location','best');
%D = f1.PaperPosition; % Returns 1x4 vector [left right width height]
%f1.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
% print -f1 -dpsc2 sets_drift1.eps;
print(gcf,'-depsc2','-loose','sets_drift1_xi2.eps');
%print('-depsc2', 'sets_drift1.eps')

f2=figure;
set(gcf,'paperpositionmode','auto')
plot(mu_z-k.*zfun,'Linewidth',1.6);hold on;
plot(mu_z-k.*zfun+hhfun2_infty(v_infty),'Linewidth',1.6);
plot(mu_z-k.*zfun+hhfun2_20(v_20),'Linewidth',1.6);
xlabel('z',FS,14);
legend2 = legend('original model','\theta=\infty, \xi=0.01','\theta=27.88, \xi=0.01');
set(legend2,FS,14, 'Location','best');
print(gcf, '-depsc2', '-loose','sets_drift2_xi2.eps')

% f3=figure;
% plot(-1.*rfun1_infty(v_infty),'Linewidth',1.6);hold on;
% %plot(-hfun1_infty(v_infty)+rfun1_infty(v_infty),'Linewidth',1.6);
% plot(-1.*rfun1_20(v_20),'Linewidth',1.6);
% plot(-hfun1_20(v_20)+rfun1_20(v_20),'Linewidth',1.6);
% xlabel('z',FS,14);
% legend1 = legend('-r^*(1),  \theta=\infty',...
% 	'-r^*(1), \theta=20', 'r^*(1)-h^*(1), \theta=20');
% set(legend1,FS,14, 'Location','best');
% print('-depsc2', '-loose','sets_r1.eps')

% f4=figure;
% plot(-rfun2_infty(v_infty),'Linewidth',1.6);hold on;
% %plot(-hfun2_infty(v_infty)+rfun2_infty(v_infty),'Linewidth',1.6);
% plot(-rfun2_20(v_20),'Linewidth',1.6);
% plot(-hfun2_20(v_20)+rfun2_20(v_20),'Linewidth',1.6);
% xlabel('z',FS,14);
% legend1 = legend('-r^*(2),  \theta=\infty',...
% 	'-r^*(2), \theta=20', 'r^*(2)-h^*(2), \theta=20');
% set(legend1,FS,14, 'Location','best');%print('-depsc2','-loose', 'sets_r2.eps')

figure;plot(-hfun1_infty(v_infty), 'Linewidth',1.6);
hold on;;plot(-hfun2_infty(v_infty), 'Linewidth',1.6);
plot(-hfun1_20(v_20), 'Linewidth',1.6);
plot(-hfun2_20(v_20), 'Linewidth',1.6);
xlabel('z',FS,14);
legend1 = legend('-h^*(1),  \theta=\infty, \xi=0.01',...
	'-h^*(2),  \theta=\infty, \xi=0.01','-h^*(1), \theta=27.88, \xi=0.01','-h^*(2), \theta=27.88, \xi=0.01');
set(legend1,FS,14, 'Location','best');
print('-depsc2', 'sets_h_xi2.eps')


xgrid = linspace(-1,1,5000);
onefun = chebfun(@(x) 1, d);
mu = mu_z-k.*zfun(xgrid);
sigma = norm(sigma_z).*onefun(xgrid);

dxc = [xgrid(2)-xgrid(1) 0.5*(xgrid(3:end)-xgrid(1:end-2)) xgrid(end)-xgrid(end-1)];
dsigma = [sigma(2)-sigma(1) 0.5*(sigma(3:end)-sigma(1:end-2)) sigma(end)-sigma(end-1)]./dxc;

dlogq = (mu - dsigma.*sigma)*2 ./ (sigma.^2) .* dxc;
logq = [0 0.5*cumsum(dlogq(1:end-1)+dlogq(2:end))];
logq = logq - max(logq) + 1; % normalization
% q = exp(logq)/sum(exp(logq).*dxc - 0.5*exp(logq(1))*dxc(1) - 0.5*exp(logq(end))*dxc(end));
q = exp(logq)/sum(0.5*sum(exp(logq(1:end-1)).*dxc(1:end-1) + exp(logq(2:end)).*dxc(2:end)));

% figure; 
% plot(xgrid,q, 'Linewidth',1.6);
% xlabel('z',FS,14);
% %print('-depsc2','-loose', 'sets_density.eps')

qfun =  chebfun.spline(xgrid,q);	

d1 = [ min(hfun1_infty(v_infty)) max(hfun1_infty(v_infty))];
d2 = [ min(hfun2_infty(v_infty)) max(hfun2_infty(v_infty))];
d3 = [ min(hfun1_20(v_20)) max(hfun1_20(v_20))-0.003];
d4 = [ min(hfun2_20(v_20))  max(hfun2_20(v_20))-0.003];
xx1 = linspace(d1(1),d1(2),10000);xx2 = linspace(d2(1),d2(2),10000);
xx3 = linspace(d3(1),d3(2),10000);xx4 = linspace(d4(1),d4(2),10000);

%zh1 = inv(hfun1_infty(v_infty),'splitting','on');%,'EPS', 1e-4);
%zh2 = inv(hfun2_infty(v_infty),'splitting','on','EPS', 1e-4);
zh1 = inv(chebfun(@(x) feval(hfun1_infty(v_infty),x), d, 'eps', 1e-6, 'splitting','on'),'eps', 1e-3,'splitting','on');
disp('zh1 computed');
zh2 = inv(chebfun(@(x) feval(hfun2_infty(v_infty),x), d, 'eps', 1e-6, 'splitting','on'),'eps', 1e-3,'splitting','on');
disp('zh2 computed');
zh3 = inv(chebfun(@(x) feval(hfun1_20(v_20),x), d, 'eps', 1e-6, 'splitting','on'),'eps', 1e-3,'splitting','on');
disp('zh3 computed');
%zh3 = inv(hfun1_20(v_20),'splitting','on');%,'EPS');%, 1e-3);
% zh4 = inv(hfun2_20(v_20),'splitting','on','EPS', 1e-4);
% zh4 = inv(chebfun.spline(xx4,feval(hfun2_20(v_20),xx4)),'splitting','on','EPS', 1e-4);
 zh4 = inv(chebfun(@(x) feval(hfun2_20(v_20),x), d, 'eps', 1e-6, 'splitting','on'),'eps', 1e-3,'splitting','on');
disp('zh4 computed');

g_infty1 = feval(qfun,zh1(xx1)).*abs(feval(diff(zh1,1),xx1));
g_infty2 = feval(qfun,zh2(xx2)).*abs(feval(diff(zh2,1),xx2));
g_20_1 = feval(qfun,zh3(xx3)).*abs(feval(diff(zh3,1),xx3));
g_20_2 = feval(qfun,zh4(xx4)).*abs(feval(diff(zh4,1),xx4));


% figure; plot(xx1,g_infty1,'.','Linewidth',1.6);hold on; 
% plot(xx2,g_infty2,'.','Linewidth',1.6);
% plot(xx3,g_20_1,'.','Linewidth',1.6);
% plot(xx4,g_20_2,'.','Linewidth',1.6);
% legend1 = legend('p(h_1),  \theta=\infty',...
% 	'p(h_2),  \theta=\infty','p(h_1), \theta=38,\xi=0.01','p(h_2), \theta=38,\xi=0.01');
% set(legend1,FS,14, 'Location','best');
% xlabel('h',FS,14);
% print('-depsc2', 'sets_h_density_ksi1.eps')

clear alpha;
figure; 
h2=area(xx3,g_20_1);hold on;
alpha(0.5);h1 = area(xx1,g_infty1);
legend1 = legend('p(h_1), \theta=27.88, \xi=0.01','p(h_1), \theta=\infty, \xi=0.01');
set(legend1,FS,14, 'Location','best');
xlabel('h',FS,14);
h1.FaceColor = [0 0.25 0.25];
%h1.EdgeColor = [0 0.25 0.25];
h1.EdgeColor = 'none';
h2.FaceColor = [0 0.5 0.5];
%h2.EdgeColor = [0 0.5 0.5];
h2.EdgeColor = 'none';
print('-depsc2','-loose', 'sets_h1_density_xi2.eps')

figure; 
h2=area(xx4,g_20_2);hold on;h1 = area(xx2,g_infty2);
legend1 = legend('p(h_2), \theta=27.88, \xi=0.01','p(h_2), \theta=\infty, \xi=0.01');
set(legend1,FS,14, 'Location','best');
xlabel('h',FS,14);
h1.FaceColor = [0 0.25 0.25];
%h1.EdgeColor = [0 0.25 0.25];
h1.EdgeColor = 'none';
h2.FaceColor = [0 0.5 0.5];
%h2.EdgeColor = [0 0.5 0.5];
h2.EdgeColor = 'none';
print('-depsc2','-loose', 'sets_h2_density_xi2.eps')

figure;
set(gcf,'paperpositionmode','auto')
plot(diff(v_20,1),'Linewidth',1.6);hold on;
plot(diff(v_infty,1),'Linewidth',1.6);
xlabel('z',FS,14); ylabel('v^{\prime}(z)',FS,14);
legend2 = legend('\theta=\infty, \xi=0.01','\theta=27.88, \xi=0.01');
set(legend2,FS,14, 'Location','best');
print(gcf, '-depsc2', '-loose','sets_value_diff_xi2.eps')
