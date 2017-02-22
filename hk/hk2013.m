load baselinesolution.mat;

chebfunpref.setDefaults('chebfuneps',1e-6);
domain_x = [0.00001, 1];
mucheb = polyfit((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),miux(2:totallength-x),57,domain(domain_x));
sigmacheb = polyfit((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),sigmax(2:totallength-x),57,domain(domain_x));
price_dividend = polyfit((F(2:totallength-x)-y(2:totallength-x))./F(2:totallength-x),F(2:totallength-x),57,domain(domain_x));

xfun = chebfun(@(x) x,domain_x);
w =  1+l-rho*(1-xfun).*price_dividend;
xi =  -rho*(diff(price_dividend,1).*(1-xfun)-price_dividend)./w;
clearvars -except domain_x mucheb sigmacheb xi g sigma rho gama

betacheb = chebfun(mucheb.*xi+g -sigma^2/2, domain_x, 'splitting', 'on'); 
alphacheb = chebfun(sigmacheb.*xi+sigma, domain_x,'splitting', 'on'); 

betacheb_sdf =  -rho - gama.*betacheb; 
alphacheb_sdf = -gama*alphacheb; 

l.boundary_type = 'natural';
r.boundary_type = 'natural';
m_name ='consumption';
sh_type1 ='shock-exposure elasticity, state space: experts wealth share';
sh_type2 ='shock-price elasticity, state space: experts wealth share';

save('hk2013.mat')