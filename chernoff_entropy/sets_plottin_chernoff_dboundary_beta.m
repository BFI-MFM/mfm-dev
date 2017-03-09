clear all; close all;
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';

set_path_tensor

load hl_theta_beta_960_80.mat
hfun1_80 = hfun1;
hfun2_80 = hfun2;
rfun1_80 = rfun1;
rfun2_80 = rfun2;
valfun_80 = valfun;
f_80 = f;
z_sol_80 = z_sol;
load  hl_theta_beta_960_120.mat
hfun1_120 = hfun1;
hfun2_120 = hfun2;
rfun1_120 = rfun1;
rfun2_120 = rfun2;
valfun_120 = valfun;
f_120 = f;
z_sol_120 = z_sol;
zzz =linspace(-0.5,0.5, 500);
zzz1=linspace(-0.1, -0.01,100);
zzz2=linspace(0, 0.5,500);
figure

plot(zz1,ff1,'Linewidth',2.6); hold on;
plot(zz2,ff2,'Linewidth',2.6); 
plot(feval(z_sol1,zzz1),feval(f1,feval(z_sol1,zzz1)),'k','Linewidth',3.6);
plot(feval(z_sol2,zzz2),feval(f2,feval(z_sol2,zzz2)),'k','Linewidth',3.6);
xlabel('\kappa', FS, 14); ylabel('\beta(\kappa)',FS, 14); 
print(gcf, '-depsc2', '-loose','sets_beta_960_3_boundary.eps')

% FS = 'fontsize';
% set(gcf,'paperpositionmode','auto')
% f1= figure;
% plot(zzz,sigma_z(1).*feval(hfun1_80,zzz)+sigma_z(2).*feval(hfun2_80,zzz),'Linewidth',1.6');hold on;
% plot(zzz,sigma_z(1).*feval(hfun1_120,zzz)+sigma_z(2).*feval(hfun2_120,zzz),'Linewidth',1.6');
% plot(zzz,sigma_z(1).*feval(hfun1_120,zzz)+sigma_z(2).*feval(hfun2_120,zzz),'Linewidth',1.6');
% plot(zzz,sigma_z(1).*feval(rfun1_80,zzz)+sigma_z(2).*feval(rfun2_80,zzz),'Linewidth',1.6');hold on;
% plot(zzz,sigma_z(1).*feval(rfun1_120,zzz)+sigma_z(2).*feval(rfun2_120,zzz),'Linewidth',1.6');
% plot(zzz,sigma_z(1).*feval(rfun1_120,zzz)+sigma_z(2).*feval(rfun2_120,zzz),'Linewidth',1.6');
% xlabel('z',FS,14);
% leg = legend('$\sigma_z\cdot H, T=80$','$\sigma_z\cdot H, T=120 $','$\sigma_z\cdot H, T=120$',...
% 		'$\sigma_z\cdot R, T=80$','$\sigma_z\cdot R, T=120 $','$\sigma_z\cdot R, T=120$');
% set(leg,'Interpreter','latex');
% set(leg,FS,14, 'Location','best');
% %D = f1.PaperPosition; % Returns 1x4 vector [left right width height]
% %f1.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
% % print -f1 -dpsc2 sets_drift1.eps;
% print(gcf,'-depsc2','-loose','sets_alpha_z_960_4_h.eps');
% %print('-depsc2', 'sets_drift1.eps')

figure; plot(zzz,feval(valfun_80,zzz),  'Linewidth',1.6);hold on; 
plot(zzz,feval(valfun_120,zzz),'Linewidth',1.6); 
xlabel('z', FS, 14); ylabel(' v(z)',FS, 14); 
leg = legend('$T=24$','$T=40 $');
set(leg,'Interpreter','latex');
set(leg,FS,14, 'Location','best');
print(gcf, '-depsc2', '-loose','sets_beta_960_3_val.eps')

figure; plot(zzz,feval(diff(valfun_80,1),zzz),  'Linewidth',1.6);hold on; 
plot(zzz,feval(diff(valfun_120,1),zzz),'Linewidth',1.6); 
xlabel('z', FS, 14); ylabel('d v/ dz',FS, 14); 
leg = legend('$T=24$','$T=40 $');
set(leg,'Interpreter','latex');
set(leg,FS,14, 'Location','best');
print(gcf, '-depsc2', '-loose','sets_beta_960_3_val_derivative.eps')

figure; plot(zzz,feval(diff(valfun_80,2),zzz), 'Linewidth',1.6);hold on; 
plot(zzz,feval(diff(valfun_120,2),zzz),'Linewidth',1.6); 
xlabel('z', FS, 14); ylabel('d2 v/ dz^2',FS, 14); 
leg = legend('$T=80$','$T=120 $');
set(leg,'Interpreter','latex');
set(leg,FS,14, 'Location','best');
print(gcf, '-depsc2', '-loose','sets_beta_960_3_val_derivative_second.eps')


FS = 'fontsize';
set(gcf,'paperpositionmode','auto')
figure;
plot(zzz,feval(-hfun1_80,zzz),'Linewidth',1.6);hold on;plot(zzz,feval(-hfun1_120,zzz),'Linewidth',1.6);
plot(zzz,feval(-rfun1_80,zzz),'Linewidth',1.6);%hold on;plot(zzz,feval(rfun1_120,zzz),'Linewidth',1.6');plot(zzz,feval(rfun1_120,zzz),'Linewidth',1.6');
xlabel('z',FS,14);
leg = legend('$-h1, T=24$','$-h1, T=40 $',...
		'$-r1, T=40$'  );
set(leg,'Interpreter','latex');
set(leg,FS,14, 'Location','best');
print(gcf,'-depsc2','-loose','sets_beta_960_3_h1.eps');

figure;
plot(zzz,feval(-hfun2_80,zzz),'Linewidth',1.6);hold on;plot(zzz,feval(-hfun2_120,zzz),'Linewidth',1.6);
plot(zzz,feval(-rfun2_80,zzz),'Linewidth',1.6);
%plot(zzz,feval(rfun2_120,zzz),'Linewidth',1.6');plot(zzz,feval(rfun2_120,zzz),'Linewidth',1.6');
xlabel('z',FS,14);
leg = legend('$-h2, T=24$','$-h2, T=40 $',...
		'$-r2, T=80$');
set(leg,'Interpreter','latex');
set(leg,FS,14, 'Location','best');
%D = f1.PaperPosition; % Returns 1x4 vector [left right width height]
%f1.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
% print -f1 -dpsc2 sets_drift1.eps;
print(gcf,'-depsc2','-loose','sets_beta_960_3_h2.eps');
%print('-depsc2', 'sets_drift1.eps')


FS = 'fontsize';
set(gcf,'paperpositionmode','auto')
figure;plot(zzz,feval(hfun1_80,zzz).*sigma_y(1)+feval(hfun2_80,zzz).*sigma_y(2)...
	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6);hold on;
plot(zzz,feval(hfun1_120,zzz).*sigma_y(1)+feval(hfun2_120,zzz).*sigma_y(2)...
	+mu_y_hat+beta_hat.*feval(zfun,zzz),'Linewidth',1.6);
plot(zzz,feval(onefun,zzz).*(mu_y_hat+beta_hat.*feval(zfun,zzz)),'Linewidth',1.6);
xlabel('z',FS,14);
lab=ylabel('$\sigma_y \cdot H +\hat{\alpha_y}+\hat{\beta} z$',FS,14);
 set(lab,'Interpreter','latex');
legend1 = legend('T=24','T=40','baseline');
set(legend1,FS,14, 'Location','best');
print(gcf,'-depsc2','-loose','sets_beta_960_3_drift1.eps');


figure; plot(zzz,feval(f_80.*zfun+ mu_y_hat.*onefun,zzz), 'Linewidth',1.6); hold on;
 plot(zzz,feval(f_120.*zfun+ mu_y_hat.*onefun,zzz), 'Linewidth',1.6); 
 lab=ylabel('$\hat{\alpha_y}+\beta z$',FS,14);
 set(lab,'Interpreter','latex');
legend1 = legend('T=24','T=40');
set(legend1,FS,14, 'Location','best');
xlabel('z', FS, 14);
print(gcf, '-depsc2', '-loose','sets_beta_960_3_drift3.eps')


FS = 'fontsize';
set(gcf,'paperpositionmode','auto')
figure;plot(zzz,feval(hfun1_80,zzz).*sigma_z(1)+feval(hfun2_80,zzz).*sigma_z(2)...
	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6);hold on;
plot(zzz,feval(hfun1_120,zzz).*sigma_z(1)+feval(hfun2_120,zzz).*sigma_z(2)...
	+mu_z_hat-kappa_hat.*feval(zfun,zzz),'Linewidth',1.6);
plot(zzz,feval(onefun,zzz).*(mu_z_hat-kappa_hat.*feval(zfun,zzz)),'Linewidth',1.6);
xlabel('z',FS,14);lab=ylabel('$\sigma_z \cdot H+\hat{\alpha_z}-\hat{\kappa} z$',FS,14);
 set(lab,'Interpreter','latex');
legend1 = legend('T=24','T=40','baseline');
set(legend1,FS,14, 'Location','best');
print(gcf,'-depsc2','-loose','sets_beta_960_3_drift2.eps');

figure; plot(zzz,feval(-z_sol_80.*zfun,zzz), LW,1.6);hold on;
	 plot(zzz,feval(-z_sol_120.*zfun,zzz), LW,1.6);hold on;
	 plot(zzz,feval(-kappa_hat.*zfun, zzz), LW,1.6);
	  xlabel('z', FS, 14); 
	leg=legend('$-\kappa(z) z, T=24$', '$-\kappa(z) z, T=40$','$-\hat{\kappa} z$','Location', 'best'); 
	set(leg,FS, 14); set(leg,'Interpreter','latex');
print(gcf, '-depsc2', '-loose','sets_beta_960_3_kappa_drift.eps')

figure; plot(zzz,feval(f_80,zzz), '.',LW,1.6);hold on;plot(zzz,feval(f_120,zzz), '.',LW,1.6);
	 plot(zzz,feval(beta_hat.*onefun,zzz), LW,1.6);
	  xlabel('z', FS, 14); leg=legend('$\beta(z), T=80$','$\beta(z), T=120$','$\hat{\beta}$', 'Location', 'best'); 
	  set(leg,FS, 14); set(leg,'Interpreter','latex');
	print(gcf, '-depsc2', '-loose','sets_beta_960_3_beta(z).eps')
