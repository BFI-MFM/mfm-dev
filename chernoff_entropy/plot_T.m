%colormap copper
set(0,'DefaultAxesColorOrder',jet(11))'  
 figure; plot(mu_z,log(2)./rr(1,:),'Linewidth',1.6);
 hold on;
plot(mu_z,log(2)./rr(50,:),'Linewidth',1.6);
 plot(mu_z,log(2)./rr(100,:),'Linewidth',1.6);
plot(mu_z,log(2)./rr(150,:),'Linewidth',1.6);
plot(mu_z,log(2)./rr(200,:),'Linewidth',1.6);
 plot(mu_z,log(2)./rr(250,:),'Linewidth',1.6);
 plot(mu_z,log(2)./rr(300,:),'Linewidth',1.6);
 plot(mu_z,log(2)./rr(350,:),'Linewidth',1.6);
 plot(mu_z,log(2)./rr(400,:),'Linewidth',1.6);
 plot(mu_z,log(2)./rr(450,:),'Linewidth',1.6);
plot(mu_z,log(2)./rr(500,:),'Linewidth',1.6);
plot(mu_z,ones(1,length(mu_z)).*160,'-k', 'Linewidth',3);
       xlabel('\mu_z',FS,14); ylabel('T',FS,14);
 legend('\kappa=1.0000e-04','\kappa=0.0081','\kappa=0.0161','\kappa=0.0241',...
'\kappa=0.0321','\kappa= 0.0401', '\kappa=0.0481','\kappa=0.0561','\kappa=0.0641', '\kappa=0.0722', '\kappa=0.08');
 set(legend,FS,14, 'Location','best');

 
 
 
 
