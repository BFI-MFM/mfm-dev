clear all; close all;
tic;
c = parcluster('local');
c.NumWorkers = 16;
par_obj= parpool(c, c.NumWorkers);
if par_obj.NumWorkers>0
    display(['working in PARALLEL mode with ',num2str(matlabpool('size')),' workers']);
else
    profile clear
    profile on
    display('working in SERIAL mode');
end
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
set_path_tensor
load beta_chernoff_set_480.mat;

mu_y_hat = 0.386;
mu_z_hat =  0;
beta_hat = 1;
kappa_hat = 0.019;
sigma_y = [0.488 0]';
sigma_z = [0.013 0.028]';
delta = 0.002;

sigma = [sigma_y'; sigma_z'];
sigma_inv = inv(sigma);
ns = 60;
chernoff_s = linspace(0.1,0.9,ns);
d1 = [zz1(1)  zz1(end)];
nrank = 24;
f1 = polyfit(zz1, ff1',nrank,domain(d1));
kappa_cheb1 = chebfun(@(z) z, d1);
d2 = [zz2(1)  zz2(end)];
f2 = polyfit(zz2, ff2',nrank,domain(d2));
iii = find(ff2==min(ff2));
k2_min = zz2(iii);
n_theta = 50;
%n_theta = 1;
%theta_space = 0.1178;
theta_space = 0.18;
theta_space = linspace(5,0.1,n_theta);
ss_sigma  = sigma*sigma';
dom_z1 = [-2 0];dom_z2 = [0 5];
xx1 = linspace(-1, -0.0001,12);
xx2 = linspace(0.0001, 1,12);
nz = 1000;
zzz1 = linspace(-2, -0.0001,nz);
dz1 = (zzz1(end)-zzz1(1))/nz;
zzz2 = linspace(0.0001,5,nz);
dz2 = (zzz2(end)-zzz2(1))/nz;

onefun1 = chebfun(@(z) 1, dom_z1, 'eps', 1e-6, 'vectorize','splitting', 'on');
zfun1 = chebfun(@(z) z, dom_z1, 'eps', 1e-6, 'vectorize','splitting', 'on');
kappa_cheb1 = chebfun(@(z) z, d1, 'eps', 1e-6, 'vectorize','splitting', 'on');
onefun2 = chebfun(@(z) 1, dom_z2, 'eps', 1e-6, 'vectorize','splitting', 'on');
zfun2 = chebfun(@(z) z, dom_z2, 'eps', 1e-6, 'vectorize','splitting', 'on');
kappa_cheb2 = chebfun(@(z) z, d2, 'eps', 1e-6, 'vectorize','splitting', 'on');
valfun1 = chebfun(@(z) 0.1.*(0.3.*z-2).^3, dom_z1, 'eps', 1e-6, 'vectorize');
valfun2 = chebfun(@(z) 0.001.*(0.3.*z-2).^3, dom_z2, 'eps', 1e-6, 'vectorize');
for kk=1:n_theta
    valfun1 = chebfun(@(z) 0.1.*(0.3.*z-2).^3, dom_z1, 'eps', 1e-6, 'vectorize');
    valfun2 = chebfun(@(z) 0.001.*(0.3.*z-2).^3, dom_z2, 'eps', 1e-6, 'vectorize');
    disp('computing theta');
    %kk=1;
    theta = theta_space(kk);
    disp(theta);
    %valfun1 = chebfun(@(z) 0.1.*(0.3.*z-2).^3, dom_z1, 'eps', 1e-6, 'vectorize');
    %valfun2 = chebfun(@(z) 0.01.*(0.3.*z-2).^3, dom_z2, 'eps', 1e-6, 'vectorize');
    %valfun = chebfun({@(z) b_minus.*z,@(z)b_plus.*z }, [dom_z(1)  0 dom_z(end)], ...
    %        'eps', 1e-4, 'vectorize');
    %valfun= w{3000};
    w1{1} = valfun1;v1 = valfun1;
    w2{1} = valfun2;v2 = valfun2;
    tol_val1 = 1; i=1;
    
    %-delta.*diff(valfun,1)+0.01.*beta -kappa;
    kappa1 = roots((-delta-kappa_cheb1).*diff(f1,1) + f1);
    if isempty(kappa1)
        kappa1 = min(zz1);
    end
    beta1 = f1(kappa1);
  
    val_dif1 =  0.01*(feval(diff(f1,1),kappa1));
    
    kappa2 = roots((-delta-kappa_cheb2).*diff(f2,1) + f2);
    beta2 = f2(kappa2);
    val_dif2 =  0.01*(feval(diff(f2,1),kappa2));
    
    %val_dif_central1 = 0.4;
    %val_dif_central2 = 0.5;
    zero_diff = (val_dif1+val_dif2)/2;
    val_zero_diff = 100;
    while (tol_val1 > 0.01)
        while (abs(val_zero_diff) > 0.01)
            val_dif_central1 = zero_diff;
            val_dif_central2 = zero_diff;
            i=i+1;
            %FOC conditon
            kappa_fun1 = @(z)  roots(0.01.*diff(f1,1)-feval(diff(valfun1,1), z));
            kappa_fun2 = @(z)  roots(0.01.*diff(f2,1)-feval(diff(valfun2,1), z));
            parfor kk = 1:nz
                if  isempty(kappa_fun1(zzz1(kk))) 
                    kk1_z(kk) = min(zz1);
                else
                    kk1_z(kk) = kappa_fun1(zzz1(kk));
                end
                if  isempty(kappa_fun2(zzz2(kk))) 
                    kk2_z(kk) = max(zz2);
                else
                    kk2_z(kk) = kappa_fun2(zzz2(kk));
                end

            end
            kappa_fun1_z = polyfit(zzz1,kk1_z,nrank,domain(dom_z1));
            kappa_fun2_z = polyfit(zzz2,kk2_z,nrank,domain(dom_z2));
            kz1 = chebfun(@(z) feval(kappa_fun1_z,z) ,dom_z1, ...
                'eps', 1e-4, 'vectorize','splitting', 'on');
            fz1 = chebfun(@(z) feval(f1,feval(kappa_fun1_z,z)) ,dom_z1, ...
                'eps', 1e-4, 'vectorize','splitting', 'on');

            kz2 = chebfun( @(z) feval(kappa_fun2_z,z) ,dom_z2, ...
                'eps', 1e-4, 'vectorize','splitting', 'on');
            fz2 = chebfun(@(z) feval(f2,feval(kappa_fun2_z,z)) ,dom_z2, ...
                'eps', 1e-4, 'vectorize','splitting', 'on');
            N1 = chebop(dom_z1);
            
            N1.lbc =  @(v) diff(v,1)-val_dif1;
            N1.rbc =  @(v) diff(v,1)-val_dif_central1;
            
            theta_term1 = @(v)  (0.01.*(0.01*ss_sigma(1,1)+diff(v,1).*ss_sigma(1,2))+...
                diff(v,1).*(0.01.*ss_sigma(2,1)+diff(v,1).*ss_sigma(2,2)))...
                *1./(2.*theta);
            N1.op =  @(v) -delta.*v + 0.01.*(mu_y.*onefun1+fz1.*zfun1)+...
                diff(v,1).*(mu_z.*onefun1-kz1.*zfun1)+...
                1/2.*norm(sigma_z).^2.*diff(v,2) - theta_term1(v);
            
            v1 =   N1\0;
            
            N2 = chebop(dom_z2);
            N2.rbc =  @(v) diff(v,1)-val_dif2;
            N2.lbc =  @(v) diff(v,1)-val_dif_central2;
            
            theta_term2 = @(v)  (0.01.*(0.01*ss_sigma(1,1)+diff(v,1).*ss_sigma(1,2))+...
                diff(v,1).*(0.01.*ss_sigma(2,1)+diff(v,1).*ss_sigma(2,2)))...
                *1./(2.*theta);
            
            N2.op =  @(v) -delta.*v + 0.01.*(mu_y*onefun2+fz2.*zfun2)+...
                diff(v,1).*(mu_z*onefun2-kz2.*zfun2)+...
                1/2.*norm(sigma_z).^2.*diff(v,2)- theta_term2(v);
            %delta_f2 = N2(w2{i-1});
            v2 =   N2\0;
            delta_f1 = v1-w1{i-1};
            delta_f2 = v2-w2{i-1};
            w2{i} = v2;
            w1{i} = v1;
            
            valfun1 = w1{i};
            tol_val1 =(abs(max(delta_f1))+abs(min(delta_f1 )))./abs(max(v1));
            disp(i); disp(tol_val1);
            valfun2 = w2{i};
            tol_val2 =(abs(max(delta_f2 ))+abs(min(delta_f2 )))./abs(max(v2));
            disp(i); disp(tol_val2);
            val_zero_diff = feval(valfun1,0)-feval(valfun2,0);
            disp(strcat('valfun level difference:', num2str(val_zero_diff)));
            
            if val_zero_diff > 0
                zero_diff = zero_diff - 0.25.*val_zero_diff;
            else
                zero_diff = zero_diff - 0.25.*val_zero_diff;
            end
            disp(strcat('derivative at zero:', num2str(zero_diff)));
        end
    end
    %disp(i); disp((abs(max(delta_f1 ))+abs(min(delta_f1 )))./abs(max(v1)));
    toc;
    
    zplot1= linspace(-1,0,500);zplot2= linspace(0,1,500);
    % FS = 'fontsize';
    % set(gcf,'paperpositionmode','auto')
    % figure;
    % plot(zplot1,feval(valfun1, zplot1),'Linewidth',1.6);hold on;plot(zplot2,feval(valfun2,zplot2),'Linewidth',1.6);
    % xlabel('z',FS,14);
    % ylabel('v(z)',FS,14);
    % leg = legend('$z<0$','$z>0$'  );
    % set(leg,'Interpreter','latex');
    % set(leg,FS,14, 'Location','best');
    % print(gcf,'-depsc2','-loose','case3_valfun.eps');
    
    
    % FS = 'fontsize';
    % set(gcf,'paperpositionmode','auto')
    % figure;
    % plot(zplot1,feval(diff(valfun1,1), zplot1),'Linewidth',1.6);hold on;plot(zplot2,feval(diff(valfun2,1),zplot2),'Linewidth',1.6);
    % xlabel('z',FS,14);
    % ylabel('v^\prime(z)',FS,14);
    % leg = legend('$z<0$','$z>0$'  );
    % set(leg,'Interpreter','latex');
    % set(leg,FS,14, 'Location','best');
    % print(gcf,'-depsc2','-loose','case3_valfun_diff.eps');
    
    dom_z = [dom_z1(1) 0 dom_z2(end)];
    valfun = chebfun({@(z) valfun1(z),@(z) valfun2(z)}, [dom_z1(1) 0 dom_z2(end)]);
    kappa_z = chebfun({@(z) kz1(z),@(z) kz2(z)}, [dom_z1(1) 0 dom_z2(end)], 'splitting','on');
    beta_z = chebfun({@(z) fz1(z),@(z) fz2(z)}, [dom_z1(1) 0 dom_z2(end)], 'splitting','on');
    zfun = chebfun(@(z) z, dom_z);
    
    rfun1 = chebfun(@(z) sigma_inv(1,1).*(mu_y-mu_y_hat)+sigma_inv(1,2).*(mu_z-mu_z_hat)+...
        sigma_inv(1,1).*(beta_z(z)-beta_hat).*zfun(z)+sigma_inv(1,2).*(kappa_hat-kappa_z(z)).*zfun(z),...
        [dom_z1(1) 0 dom_z2(end)], 'splitting','on');
    rfun2 = chebfun(@(z) sigma_inv(2,1).*(mu_y-mu_y_hat)+sigma_inv(2,2).*(mu_z-mu_z_hat)+...
        sigma_inv(2,1).*(beta_z(z)-beta_hat).*zfun(z)+sigma_inv(2,2).*(kappa_hat-kappa_z(z)).*zfun(z),...
        [dom_z1(1) 0 dom_z2(end)], 'splitting','on');
    hfun1 =  chebfun(@(z) rfun1(z) - 1./theta.*(sigma(1,1).*0.01+sigma(2,1).*feval(diff(valfun,1),z)),...
        [dom_z1(1) 0 dom_z2(end)], 'splitting','on');
    hfun2 =  chebfun(@(z) rfun2(z) - 1./theta.*(sigma(1,2).*0.01+sigma(2,2).*feval(diff(valfun,1),z)),...
        [dom_z1(1) 0 dom_z2(end)], 'splitting','on');
    %figure;plot(-hfun1,LW,1.6);hold on;plot(-hfun2,'r',LW,1.6); xlabel('z', FS, 14); ylabel('-h^*(z)',FS, 14);
    %h=legend('-h_1(z)','-h_2(z)','Location', 'best'); set(h,FS, 14);
    %file_save = ['sets_density-theta-xi-0-01-' num2str(kk) '.mat'];
    %save(file_save);
    
    parfor jj = 1:ns
        N_c = chebop(dom_z);
        s = chernoff_s(jj);
        tic;
        N_c.op =    @(e) -0.5.*s.*(1-s).*e.*(hfun1.^2+hfun2.^2) +...
            (s.*(sigma_z(1).*hfun1+sigma_z(2).*hfun2)+ mu_z_hat-kappa_hat.*zfun).*diff(e,1)+...
            1/2.*norm(sigma_z).^2.*diff(e,2);
        N_c.bc = 'neumann';
        [V,D]= eigs(N_c,100);
        e_v =  -sort((real(diag(D))),'descend');
        disp(e_v(1));
        rho1(jj,kk) = e_v(1);
        disp(jj);
        toc;
        N2_c = chebop(dom_z);
        s = chernoff_s(jj);
        tic;
        N2_c.op =   @(e) -0.5.*s.*(1-s).*e.*(rfun1.^2+rfun2.^2) +...
            (s.*(sigma_z(1).*rfun1+sigma_z(2).*rfun2) + mu_z_hat-kappa_hat.*zfun).*diff(e,1)+...
            1/2.*norm(sigma_z).^2.*diff(e,2);
        N2_c.bc = 'neumann';
        [V,D]= eigs(N2_c,100);
        e_v =  -sort((real(diag(D))),'descend');
        disp(e_v(1));
        rho2(jj,kk) = e_v(1);
        disp(jj);
        toc;
    end
    disp(strcat('theta=',num2str(theta)));
     disp(strcat('hl=',num2str(log(2)/max(rho1(:,kk)))))
    file_save = ['hl_theta_beta_480_linear-' num2str(kk) '.mat'];
    %file_save = ['hl_theta_beta_960_120.mat'];
    save(file_save);
end
file_save = ['hl_theta_beta_480-' num2str(kk) '.mat'];
save(file_save);
delete(par_obj)