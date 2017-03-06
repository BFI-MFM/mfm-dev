load baselinesolution_log.mat;

chebfunpref.setDefaults('chebfuneps',1e-4);
domain_x = [0.001, 1];
mucheb = polyfit((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),miux(2:totallength-x),17,domain(domain_x));
sigmacheb = polyfit((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),sigmax(2:totallength-x),17,domain(domain_x));
%price_dividend = polyfit((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),F(2:totallength-x),17,domain(domain_x));
price_dividend =chebfun(@(x) (1+l)/rho, domain_x); 
xfun = chebfun(@(x) x,domain_x);
w =  1+l-rho*(1-xfun).*price_dividend;
%for log function
xi=rho.*price_dividend./w;
%xi =  -rho*(diff(price_dividend,1).*(1-xfun)-price_dividend)./w;

figure(1);
plot((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),miux(2:totallength-x));
hold on;
plot((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),...
    feval(mucheb,(F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),miux(2:totallength-x)));

plot(eta,feval(sigmacheb,eta));  
xlabel('\eta')
ylabel('\sigma^\eta \eta');

clearvars -except domain_x mucheb sigmacheb xi g sigma rho gama

betacheb = chebfun(mucheb.*xi+g -sigma^2/2, domain_x, 'splitting', 'on'); 
alphacheb = chebfun(sigmacheb.*xi+sigma, domain_x,'splitting', 'on'); 

betacheb_sdf =  -rho - gama.*betacheb; 
alphacheb_sdf = -gama*alphacheb; 

l.boundary_type = 'reflecting';
r.boundary_type = 'reflecting';
m_name ='specialist consumption';
sh_type =' state space: experts wealth share';

save('hk2013_log.mat')