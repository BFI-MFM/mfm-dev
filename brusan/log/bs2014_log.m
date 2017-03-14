%vz: 03/13/2017 to transform BS model outputs into DVD inputs
load brusan_log.mat;
chebfunpref.setDefaults('chebfuneps',1e-6);
domain_x = [eta(1) eta(end)];
ee = eta(netaout);
m1 = polyfit(eta(1:netaout), mu_eta(1:netaout).*eta(1:netaout),16,domain([domain_x(1) eta(netaout)]));
m2 = polyfit(eta(1+netaout:end), mu_eta(1+netaout:end).*eta(1+netaout:end),6,domain([eta(netaout+1) domain_x(end)]));
mucheb = chebfun({@(x) feval(m1,x),@(x) feval(m2,x)},domain([domain_x(1) eta(netaout) domain_x(end)]));
s1 = polyfit(eta(1:netaout), sigma_eta(1:netaout).*eta(1:netaout),16,domain([domain_x(1) eta(netaout)]));
s2 = polyfit(eta(1+netaout:end), sigma_eta(1+netaout:end).*eta(1+netaout:end),6,domain([eta(netaout+1) domain_x(end)]));
sigmacheb = chebfun({@(x) feval(s1,x),@(x) feval(s2,x)},domain([domain_x(1) eta(netaout) domain_x(end)]));
phicheb = polyfit(eta, Phi(q),29,domain(domain_x));
r1 = polyfit(eta(1:netaout), r(1:netaout),46,domain([domain_x(1) eta(netaout)]));
r2 = polyfit(eta(1+netaout:end), r(1+netaout:end),6,domain([eta(netaout+1) domain_x(end)]));
rcheb = chebfun({@(x) feval(r1,x),@(x) feval(r2,x)},domain([domain_x(1) eta(netaout) domain_x(end)]));
q1 = polyfit(eta(1:netaout), q(1:netaout),46,domain([domain_x(1) eta(netaout)]));
q2 = polyfit(eta(1+netaout:end), q(1+netaout:end),6,domain([eta(netaout+1) domain_x(end)]));
qcheb = chebfun({@(x) feval(q1,x),@(x) feval(q2,x)},domain([domain_x(1) eta(netaout) domain_x(end)]));

figure(1);
subplot(1,2,1); hold on
plot(eta, sigma_eta.*eta);  
plot(eta,feval(sigmacheb,eta));  
xlabel('\eta')
ylabel('\sigma^\eta \eta');

subplot(1,2,2); hold on
plot(eta, mu_eta.*eta);  
plot(eta,feval(mucheb,eta)); 
xlabel('\eta')
ylabel('\mu^\eta \eta');

figure(2);
subplot(1,2,1); hold on
plot(eta, r);  
plot(eta,feval(rcheb,eta));  
xlabel('\eta')
ylabel('r');

subplot(1,2,2); hold on
plot(eta, Phi(q));  
plot(eta,feval(phicheb,eta)); 
xlabel('\eta')
ylabel('\phi');

figure(3);
hold on
plot(eta, q);  
plot(eta,feval(qcheb,eta));  
xlabel('\eta')
ylabel('q');

clearvars -except domain_x mucheb sigmacheb rcheb phicheb qcheb delta sigma rho

betacheb = chebfun(phicheb -delta -sigma^2/2, domain_x,'splitting', 'on');
alphacheb = chebfun(sigma, domain_x,'splitting', 'on');
betacheb_sdf = -rho-betacheb;
alphacheb_sdf = -alphacheb;

l.boundary_type = 'killing';
r.boundary_type = 'killing';
m_name ='agg. consumption';
sh_type ='state space: experts wealth share';

save('bs2014_log.mat')