%vz, 2016
%Brock-Mirman model with full depreciation
tic;
%setup parameters
alpha = 1./3; %capital share
beta = 0.97; %discount rate

tol_val =1000; tol_val_min = 1e-9; %accuracy
domain_k = [0 3];
kfun = chebfun(@(k) k, domain_k);
kfun0 =  0.1.*kfun.^alpha;
i = 1;
while (tol_val > tol_val_min)
    kfun= chebfun(@(k) kfun0((kfun0(k))),domain_k, 'eps', tol_val_min, 'vectorize','splitting','on');
    kprime = chebfun(@(k) kfun0(k).^(1-alpha).*kfun(k),domain_k, 'eps', tol_val_min, 'vectorize','splitting','on');
    
    dfun = chebfun(@(k) (k.^alpha+1./(alpha.*beta).*kprime(k))...
        .*alpha.*beta./(1+alpha*beta),domain_k, 'eps', tol_val_min, 'vectorize','splitting','on');
    tol_val = abs(max(dfun)-max(kfun0))./abs(max(kfun0));
    kfun0 = dfun;
    disp(i); disp(tol_val);
    i=i+1;
end
figure
plot(alpha*beta*kfun.^alpha, 'LineWidth', 2);
xlabel('$k$', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
ylabel('$k^\prime(k)$', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
title('Exact and Numerical Solutions', 'interpreter', 'latex', 'fontsize',12, 'Color', 'blue');
hold on
plot(kfun0, '--r', 'LineWidth', 3);
hold off
save('bm1972.mat')
toc;
