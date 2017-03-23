clear all; close all;
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';


load hl_theta_alpha_z_960_80.mat
hfun1_80 = hfun1;
hfun2_80 = hfun2;
rfun1_80 = rfun1;
rfun2_80 = rfun2;
load hl_theta_alpha_z_960_100.mat
hfun1_100 = hfun1;
hfun2_100 = hfun2;
rfun1_100 = rfun1;
rfun2_100 = rfun2;
load hl_theta_alpha_z_960_120.mat
hfun1_120= hfun1;
hfun2_120= hfun2;
rfun1_120= rfun1;
rfun2_120= rfun2;

load hl_theta_alpha_z_960_5-20.mat
hfun1_474= hfun1;
hfun2_474= hfun2;
rfun1_474= rfun1;
rfun2_474= rfun2;

zzz =linspace(-0.5,0.5, 500);

% FS = 'fontsize';
% set(gcf,'paperpositionmode','auto')
% f1= figure;
% plot(zzz,sigma_z(1).*feval(hfun1_80,zzz)+sigma_z(2).*feval(hfun2_80,zzz),'Linewidth',1.6');hold on;
% plot(zzz,sigma_z(1).*feval(hfun1_100,zzz)+sigma_z(2).*feval(hfun2_100,zzz),'Linewidth',1.6');
% plot(zzz,sigma_z(1).*feval(hfun1_120,zzz)+sigma_z(2).*feval(hfun2_120,zzz),'Linewidth',1.6');
% plot(zzz,sigma_z(1).*feval(rfun1_80,zzz)+sigma_z(2).*feval(rfun2_80,zzz),'Linewidth',1.6');hold on;
% plot(zzz,sigma_z(1).*feval(rfun1_100,zzz)+sigma_z(2).*feval(rfun2_100,zzz),'Linewidth',1.6');
% plot(zzz,sigma_z(1).*feval(rfun1_120,zzz)+sigma_z(2).*feval(rfun2_120,zzz),'Linewidth',1.6');
% xlabel('z',FS,14);
% leg = legend('$\sigma_z\cdot H, T=80$','$\sigma_z\cdot H, T=100 $','$\sigma_z\cdot H, T=120$',...
% 		'$\sigma_z\cdot R, T=80$','$\sigma_z\cdot R, T=100 $','$\sigma_z\cdot R, T=120$');
% set(leg,'Interpreter','latex');
% set(leg,FS,14, 'Location','best');
% %D = f1.PaperPosition; % Returns 1x4 vector [left right width height]
% %f1.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
% % print -f1 -dpsc2 sets_drift1.eps;
% print(gcf,'-depsc2','-loose','sets_alpha_z_960_4_h.eps');
% %print('-depsc2', 'sets_drift1.eps')


FS = 'fontsize';
set(gcf,'paperpositionmode','auto')
figure;
plot(zzz,feval(hfun1_80,zzz),'Linewidth',1.6');hold on;plot(zzz,feval(hfun1_100,zzz),'Linewidth',1.6');
plot(zzz,feval(hfun1_120,zzz),'Linewidth',1.6');plot(zzz,feval(hfun1_474,zzz),'Linewidth',1.6');
plot(zzz,feval(rfun1_80,zzz),'Linewidth',1.6');%hold on;plot(zzz,feval(rfun1_100,zzz),'Linewidth',1.6');plot(zzz,feval(rfun1_120,zzz),'Linewidth',1.6');
xlabel('z',FS,14);
leg = legend('$h1, T=80$','$h1, T=100 $','$h1, T=120$','$h1, T=474$',...
		'$r1, T=80$'  );
set(leg,'Interpreter','latex');
set(leg,FS,14, 'Location','best');
print(gcf,'-depsc2','-loose','sets_alpha_z_960_4_h1.eps');
figure;
plot(zzz,feval(hfun2_80,zzz),'Linewidth',1.6');hold on;plot(zzz,feval(hfun2_100,zzz),'Linewidth',1.6');
plot(zzz,feval(hfun2_120,zzz),'Linewidth',1.6');plot(zzz,feval(hfun2_474,zzz),'Linewidth',1.6');
plot(zzz,feval(rfun2_80,zzz),'Linewidth',1.6');
%plot(zzz,feval(rfun2_100,zzz),'Linewidth',1.6');plot(zzz,feval(rfun2_120,zzz),'Linewidth',1.6');
xlabel('z',FS,14);
leg = legend('$h2, T=80$','$h2, T=100 $','$h2, T=120$','$h2, T=474$',...
		'$r2, T=80$');
set(leg,'Interpreter','latex');
set(leg,FS,14, 'Location','best');
%D = f1.PaperPosition; % Returns 1x4 vector [left right width height]
%f1.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
% print -f1 -dpsc2 sets_drift1.eps;
print(gcf,'-depsc2','-loose','sets_alpha_z_960_4_h2.eps');
%print('-depsc2', 'sets_drift1.eps')


FS = 'fontsize';
set(gcf,'paperpositionmode','auto')
figure;plot(zzz,feval(hfun1_80,zzz).*sigma_y(1)+feval(hfun2_80,zzz).*sigma_y(2)...
	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6');hold on;
plot(zzz,feval(hfun1_100,zzz).*sigma_y(1)+feval(hfun2_100,zzz).*sigma_y(2)...
	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6');
plot(zzz,feval(hfun1_120,zzz).*sigma_y(1)+feval(hfun2_120,zzz).*sigma_y(2)...
	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6');
plot(zzz,feval(hfun1_474,zzz).*sigma_y(1)+feval(hfun2_474,zzz).*sigma_y(2)...
	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6');
plot(zzz,feval(onefun,zzz).*(mu_y_hat+beta_hat.*feval(zfun,zzz)),'Linewidth',1.6');
xlabel('z',FS,14);
lab=ylabel('$\sigma_y \cdot H +\hat{\alpha_y}+\hat{\beta} z$',FS,14);
 set(lab,'Interpreter','latex');
legend1 = legend('T=80','T=100','T=120','T=474','baseline');
set(legend1,FS,14, 'Location','best');
print(gcf,'-depsc2','-loose','sets_alpha_z_960_4_drift1.eps');


FS = 'fontsize';
set(gcf,'paperpositionmode','auto')
figure;plot(zzz,feval(hfun1_80,zzz).*sigma_z(1)+feval(hfun2_80,zzz).*sigma_z(2)...
	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6');hold on;
plot(zzz,feval(hfun1_100,zzz).*sigma_z(1)+feval(hfun2_100,zzz).*sigma_z(2)...
	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6');
plot(zzz,feval(hfun1_120,zzz).*sigma_z(1)+feval(hfun2_120,zzz).*sigma_z(2)...
	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6');
plot(zzz,feval(hfun1_474,zzz).*sigma_z(1)+feval(hfun2_474,zzz).*sigma_z(2)...
	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6');
plot(zzz,feval(onefun,zzz).*(mu_z_hat-kappa_hat.*feval(zfun,zzz)),'Linewidth',1.6');
xlabel('z',FS,14);lab=ylabel('$\sigma_z \cdot H+\hat{\alpha_z}-\hat{\kappa} z$',FS,14);
 set(lab,'Interpreter','latex');
legend1 = legend('T=80','T=100','T=120','T=474','baseline');
set(legend1,FS,14, 'Location','best');
print(gcf,'-depsc2','-loose','sets_alpha_z_960_4_drift2.eps');