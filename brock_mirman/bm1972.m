%vz, 2016
%Brock-Mirman model with full depreciation
tic;
%setup parameters
alpha = 1./3; %capital share
beta = 0.97; %discount rate

tol_val =1000; tol_val_k= 1e-5; tol_cheb=1e-5;%accuracy
domain_k = [0 3];
kfun = chebfun(@(k) k, domain_k);
k_pol0 =  0.1.*kfun.^alpha;
i = 1;
while (tol_val > tol_val_k)
    k_pol= chebfun(@(k) k_pol0((k_pol0(k))),domain_k, 'eps', tol_cheb, 'vectorize','splitting','on');
    kprime = chebfun(@(k) k_pol0(k).^(1-alpha).*k_pol(k),domain_k, 'eps', tol_cheb, 'vectorize','splitting','on');
    
    dfun = chebfun(@(k) (k.^alpha+1./(alpha.*beta).*kprime(k))...
        .*alpha.*beta./(1+alpha*beta),domain_k, 'eps', tol_cheb, 'vectorize','splitting','on');
    tol_val = abs(max(dfun)-max(k_pol0))./abs(max(k_pol0));
    k_pol0 = dfun;
    disp(i); disp(tol_val);
    i=i+1;
end
figure
plot(alpha*beta*kfun.^alpha, 'LineWidth', 2);
xlabel('$k$', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
ylabel('$k^\prime(k)$', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
title('Exact and Numerical Solutions', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
hold on
plot(k_pol0, '--r', 'LineWidth', 3);
legend('Exact', 'Numerical', 'Location', 'best');
hold off
xx = linspace(domain_k(1), domain_k(end),200);
figure
plot(xx,alpha*beta.*kfun(xx).^alpha-k_pol0(xx), 'LineWidth', 2);
xlabel('$k$', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
ylabel('$k^\prime(k)$', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
title('Numerical Solutions error', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
save('bm1972.mat')
toc;
