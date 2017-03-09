clear all; close all;

 % c = parcluster('local');
 % c.NumWorkers = 16;
 % par_obj= parpool(c, c.NumWorkers);
 % if par_obj.NumWorkers>0
 %     display(['working in PARALLEL mode with ',num2str(matlabpool('size')),' workers']);
 % else
 %     profile clear
 %     profile on
 %     display('working in SERIAL mode');
 % end
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
set_path_tensor
load beta_chernoff_set_960.mat;

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
f1= polyfit(zz1, ff1',nrank,domain(d1));

d2 = [zz2(1)  zz2(end)];
 f2= polyfit(zz2, ff2',nrank,domain(d2));
n_theta = 20;
n_theta = 1;
%theta_space = 0.1178;
theta_space = 100;
%theta_space = linspace(0.5,0.1,n_theta);
ss_sigma  = sigma*sigma';
dom_z1 = [-50 0];dom_z2 = [0 50];
xx1 = linspace(-1, -0.0001,12);
xx2 = linspace(0.0001, 1,12);
zzz1 = linspace(-2, -0.00001,5000);
zzz2 = linspace(0.00001,2,5000);
onefun1 = chebfun(@(z) 1, dom_z1, 'eps', 1e-6, 'vectorize','splitting', 'on');
zfun1 = chebfun(@(z) z, dom_z1, 'eps', 1e-6, 'vectorize','splitting', 'on');

onefun2 = chebfun(@(z) 1, dom_z2, 'eps', 1e-6, 'vectorize','splitting', 'on');
zfun2 = chebfun(@(z) z, dom_z2, 'eps', 1e-6, 'vectorize','splitting', 'on');

for kk=1:n_theta
    disp('computing theta');
    %kk=1;
    theta = theta_space(kk);
    disp(theta);

    %valfun = 0.05.*chebfun(@(z) z+cos(0.3.*z-1), dom_z, 'splitting', 'on','eps', 1e-6, 'vectorize');
    valfun1 = chebfun(@(z) 0.1.*(0.3.*z-2).^3, dom_z1, 'eps', 1e-6, 'vectorize');
    valfun2 = chebfun(@(z) 0.01.*(0.3.*z-2).^3, dom_z2, 'eps', 1e-6, 'vectorize');
    %valfun = chebfun({@(z) b_minus.*z,@(z)b_plus.*z }, [dom_z(1)  0 dom_z(end)], ...
    %        'eps', 1e-4, 'vectorize');
    %valfun= w{3000};
    w1{1} = valfun1;v1 = valfun1;
    w2{1} = valfun2;v2 = valfun2;
    tol_val1 = 1; i=1;

    %-delta.*diff(valfun,1)+0.01.*beta -kappa;
    while (tol_val1 > 0.01)
        i=i+1;
        %FOC conditon
        u1 = 0.01.*diff(f1,1)-feval(diff(valfun1,1), zzz1);
        temp = roots(u1);
        z_roots1 =temp(1,:);
        zzr1 = zzz1;ffr1 = ff1';
        zzr1(isnan(z_roots1)) = [];
        z_roots1(isnan(z_roots1)) = [];

        dr1 = [min(zzr1) max(zzr1)];
        if size(zzr1,2)>4
            kappa_fun1 = polyfit(zzr1,z_roots1,4,domain(dr1));
        else
            kappa_fun1 = chebfun(@(x) mean(z_roots1), domain(d1))
        end

        u2 = 0.01.*diff(f2,1)-feval(diff(valfun2,1), zzz2);
        temp = roots(u2);
        z_roots2 =temp(1,:);
        zzr2 = zzz2;ffr2 = ff2';
        zzr2(isnan(z_roots2)) = [];z_roots2(isnan(z_roots2)) = [];
        dr2 = [min(zzr2) max(zzr2)];
        if size(zzr2,2)>4
            kappa_fun2 = polyfit(zzr2,z_roots2,4,domain(dr2));
        else
            kappa_fun2 = chebfun(@(x) mean(z_roots2), domain(d2))
        end

        %if min(roots_fun) < dom_z(1)
        %   min(roots_fun) = dom_z(1);
        %end
        dom_z1 =   [kappa_fun1.domain(1) kappa_fun1.domain(2)];
        dom_z2 =   [kappa_fun2.domain(1) kappa_fun2.domain(2)];
        disp('starting z_sol');
        z_sol1 = kappa_fun1;
        z_sol2 = kappa_fun2;
        fz1.domain = [min(z_sol1) max(z_sol1)];
        fz2.domain =  [min(z_sol2) max(z_sol2)];
        fz1 = chebfun(@(z) feval(f1,feval(z_sol1,z)) , [dom_z1(1) dom_z1(2) ], ...
            'eps', 1e-5, 'vectorize');
        fz2 = chebfun(@(z) feval(f2,feval(z_sol2,z)) , [ dom_z2(1)  dom_z2(2)], ...
            'eps', 1e-5, 'vectorize');
         disp('starting PDE');
        %cross_point1(i) = feval(diff(w1{i-1},1), 0);
        %cross_point2(i) = feval(diff(w2{i-1},1), 0); 

        cross_point(i) = feval(chebfun.spline([xx1 xx2],[feval(diff(v1,1),xx1) feval(diff(v2,1),xx2)]),0);
        %valzero(i) = (feval(w1{i-1}, 0)+feval(w2{i-1},0))/2;
        valzero(i) = feval(chebfun.spline([xx1 xx2],[feval(v1,xx1) feval(v2,xx2)]),0);
        chebfun.spline([xx1 xx2(1)],[v1(xx1) v2(xx2(1))]);
        disp(strcat('value at zero left:', num2str(feval(v1,0))));
        disp(strcat('value at zero right:', num2str(feval(v2,0))));
        disp(strcat('derivative at zero left:', num2str(feval(diff(v1,1),0))));
        disp(strcat('derivative at zero right:', num2str(feval(diff(v2,1),0))));
        N1 = chebop(dom_z1);
        %N1.rbc = [valzero(i); cross_point(i)];
        %N1.rbc = [cross_point(i)]
        %N1.rbc =  @(v) diff(v,1)-cross_point;
        %N1.rbc =  @(v) [v-valzero(i);  diff(v)-cross_point];
        N1.lbc = @(v) (-delta-kappa_fun1).*diff(v,1)+0.01.*fz1;
        w1{i-1}.domain = dom_z1;   onefun1.domain = dom_z1;zfun1.domain = dom_z1;
        theta_term1 = @(v)  (0.01.*(0.01*ss_sigma(1,1)+diff(v,1).*ss_sigma(1,2))+...
            diff(v,1).*(0.01.*ss_sigma(2,1)+diff(v,1).*ss_sigma(2,2)))...
            *1./(2.*theta);
        N1.op =  @(v) -delta.*v + 0.01.*(mu_y.*onefun1+fz1.*zfun1)+diff(v,1).*(mu_z*onefun1-z_sol1.*zfun1)+...
            1/2.*norm(sigma_z).^2.*diff(v,2)- theta_term1(v);

        delta_f1 = N1(w1{i-1});
        v1 =   w1{i-1} + delta_f1;   

        N2 = chebop(dom_z2);
        N2.rbc = @(v)  (-delta-kappa_fun2).*diff(v,1)+0.01.*fz2;
        %N2.lbc =  @(v) [v-valzero(i);  diff(v)-cross_point];
        %N2.lbc =  @(v) diff(v)-cross_point;
        w2{i-1}.domain = dom_z2;   onefun2.domain = dom_z2;zfun2.domain = dom_z2;
        theta_term2 = @(v)  (0.01.*(0.01*ss_sigma(1,1)+diff(v,1).*ss_sigma(1,2))+...
            diff(v,1).*(0.01.*ss_sigma(2,1)+diff(v,1).*ss_sigma(2,2)))...
            *1./(2.*theta);
        N2.op =  @(v) -delta.*v + 0.01.*(mu_y*onefun2+fz2.*zfun2)+diff(v,1).*(mu_z*onefun2-z_sol2.*zfun2)+...
            1/2.*norm(sigma_z).^2.*diff(v,2)- theta_term2(v);
        delta_f2 = N2(w2{i-1});
        v2 =   w2{i-1} + delta_f2;
        w2{i} = v2;
        w1{i} = v1;
        %w1{i} = chebfun.spline([xx1 xx2(1)],[v1(xx1) v2(xx2(1))]);
        %w{i} = chebfun.spline(xx,v(xx));
        valfun1 = w1{i};
        tol_val1 = (abs(max(delta_f1 ))+abs(min(delta_f1 )))./abs(max(v1));
        disp(i); disp(tol_val1);
        %w2{i} = chebfun.spline([xx1(end) xx2],[v1(xx1(end)) v2(xx2)]);
        valfun2 = w2{i};
        tol_val2 = (abs(max(delta_f2 ))+abs(min(delta_f2 )))./abs(max(v2));
        disp(i); disp(tol_val2);
    end
    disp(i); disp((abs(max(delta_f1 ))+abs(min(delta_f1 )))./abs(max(v1)));

    rfun1 = sigma_inv(1,1).*(mu_y-mu_y_hat)+sigma_inv(1,2).*(mu_z-mu_z_hat)+...
        sigma_inv(1,1).*(f1-beta_hat).*zfun+sigma_inv(1,2).*(kappa_hat-z_sol).*zfun;
    rfun2 = sigma_inv(2,1).*(mu_y-mu_y_hat)+sigma_inv(2,2).*(mu_z-mu_z_hat)+...
        sigma_inv(2,1).*(f-beta_hat).*zfun+sigma_inv(2,2).*(kappa_hat-z_sol).*zfun;
    hfun1 =   rfun1 - 1./theta.*(sigma(1,1).*0.01+sigma(2,1).*diff(valfun,1));
    hfun2 =  rfun2 - 1./theta.*(sigma(1,2).*0.01+sigma(2,2).*diff(valfun,1));
    figure;plot(-hfun1,LW,1.6);hold on;plot(-hfun2,'r',LW,1.6); xlabel('z', FS, 14); ylabel('-h^*(z)',FS, 14);
    h=legend('-h_1(z)','-h_2(z)','Location', 'best'); set(h,FS, 14);
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
    file_save = ['hl_theta_beta_120_linear-' num2str(kk) '.mat'];
    %file_save = ['hl_theta_beta_960_120.mat'];
    save(file_save);
end
delete(par_obj)