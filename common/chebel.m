function [sh_el, sh_el2, pr_el, pr_el2] = chebel(model_name)
%vz: 03/13/2017 to compute DVD by using HJB PDE
%sh_el - shock elasticty of the first kind
%sh_el2 - shock elasticty of the second kind
%pr_el - price elasticty of the first kind
%pr_el2 - price elasticty of the second kind

load(model_name);
disp(model_name);

chebfunpref.setDefaults('chebfuneps',1e-4);
aa1 = chebfun(@(x) feval(1/2.*sigmacheb.^2,x), domain_x, 'splitting', 'on', 'vectorize');
bb1 = chebfun(@(x) feval(mucheb+sigmacheb.*alphacheb,x), domain_x,'splitting', 'on', 'vectorize');
cc1 = chebfun(@(x) feval(betacheb+alphacheb.^2./2,x), domain_x, 'splitting', 'on', 'vectorize');

alphacheb_val = alphacheb + alphacheb_sdf;
betacheb_val = betacheb + betacheb_sdf;
aa2 = chebfun(@(x) feval(1/2.*sigmacheb.^2,x), domain_x, 'splitting', 'on', 'vectorize');
bb2 = chebfun(@(x) feval(mucheb+sigmacheb.*alphacheb_val,x), domain_x,'splitting', 'on', 'vectorize');
cc2 = chebfun(@(x) feval(betacheb_val+alphacheb_val.^2./2,x), domain_x, 'splitting', 'on', 'vectorize');

N = chebop(domain_x);
%set left boundary conditions
if strcmp(l.boundary_type,'natural') 
	%Neumann - natural
	N.lbc =  @(v) diff(v,1);
elseif strcmp(l.boundary_type,'killing')
	%Dirichlet - killing
	N.lbc =  @(v) v-1;
elseif strcmp(l.boundary_type, 'reflecting')
	N.lbc =  @(v) diff(v,1)+v.*xi(domain_x(1));
end  

%set right boundary conditions
if strcmp(r.boundary_type,'natural') 
	%Neumann - natural
	N.rbc =  @(v) diff(v,1);
elseif strcmp(r.boundary_type,'killing')
	%Dirichlet - killing
	N.rbc =  @(v) v-1;
elseif strcmp(r.boundary_type, 'reflecting')
	N.rbc =  @(v) diff(v,1)+v.*xi(domain_x(end));
end  
  
N.op =  @(v) aa1.*diff(v,2)+bb1.*diff(v,1)+cc1.*v;
%u =  N\1;
	tic;
	u0 = chebfun(@(x) 1, domain_x);
	t = [0 0.0001 0.01 0.5 1 1.5 2 2.5 3 4 5 6 7 8 9 10 15 20 25 30 50];
	prefs = cheboppref();
	prefs.discretization = @chebcolloc2;
    prefs.plotting='on'; 
    prefs. bvpTol='1e-6'; 
	u = expm(N, t, u0, prefs);
	nx = 1000;
	xx = linspace(domain_x(1), domain_x(end),nx);
	for i=1:size(t,2)
		sh_el{i} = diff(u{i},1)./u{i}.*sigmacheb+alphacheb;%#ok<AGROW>
		sh_el_xx(i,:) = feval(sh_el{i},xx);%#ok<AGROW>
	end
	toc;
	figure;
	plot(t,sh_el_xx(:,int16(nx/25)), 'linewidth', 2);hold on;
	plot(t,sh_el_xx(:,int16(nx/2)), 'linewidth', 2); 
	plot(t,sh_el_xx(:,int16(24*nx/25)), 'linewidth', 2); 
	ylim([min(min(sh_el_xx))-0.01 max(max(sh_el_xx))+0.01]);
	xlabel('time, years'); ylabel(m_name);title(strcat(sh_type, ', type 1 elasticity, shock exposure'));
	legend(num2str(xx(int16(nx/25)), '%4.2f'),num2str(xx(int16(nx/2)), '%4.2f'),...
		num2str(xx(int16(24*nx/25)), '%4.2f'), 'Location', 'best')
%else
%	disp('state is not stochastically stable')
%end
	tic;
	u0 = chebfun(@(x) 1, domain_x);
	t = [0 0.0001 0.01 0.5 1 1.5 2 2.5 3 4 5 6 7 8 9 10 15 20 25 30 50];
	prefs = cheboppref();
	prefs.discretization = @chebcolloc2;
	u = expm(N, t, u0, prefs);
	u02 = chebfun(@(x) alphacheb(x), domain_x);
	u2 = expm(N, t, u02, prefs);
	nx = 1000;
	xx = linspace(domain_x(1), domain_x(end),nx);
	for i=1:size(t,2)
		sh_el2{i} = u2{i}./u{i};%#ok<AGROW>
		sh_el_xx2(i,:) = feval(sh_el2{i},xx);%#ok<AGROW>
	end
	toc;
	figure;
	plot(t,sh_el_xx2(:,int16(nx/25)), 'linewidth', 2);hold on;
	plot(t,sh_el_xx2(:,int16(nx/2)), 'linewidth', 2); 
	plot(t,sh_el_xx2(:,int16(24*nx/25)), 'linewidth', 2); 
	ylim([min(min(sh_el_xx2))-0.01 max(max(sh_el_xx2))+0.01]);
	xlabel('time, years'); ylabel(m_name);title(strcat(sh_type, ', type 2 elasticity, shock exposure'));
	legend(num2str(xx(int16(nx/25)), '%4.2f'),num2str(xx(int16(nx/2)), '%4.2f'),...
		num2str(xx(int16(24*nx/25)), '%4.2f'), 'Location', 'best')

N_sdf = chebop(domain_x);
%set left boundary conditions
if strcmp(l.boundary_type,'natural') 
	%Neumann - natural
	N_sdf.lbc =  @(v) diff(v,1);
elseif strcmp(l.boundary_type,'killing')
	%Dirichlet - killing
	N_sdf.lbc =  @(v) v-1;
elseif strcmp(l.boundary_type, 'reflecting')
	xi_val = xi + xi_sdf;
	N_sdf.lbc =  @(v) diff(v,1)+v.*xi_val(domain_x(1));
end  

%set right boundary conditions
if strcmp(r.boundary_type,'natural') 
	%Neumann - natural
	N_sdf.rbc =  @(v) diff(v,1);
elseif strcmp(r.boundary_type,'killing')
	%Dirichlet - killing
	N_sdf.rbc =  @(v) v-1;
elseif strcmp(r.boundary_type, 'reflecting')
	xi_val = xi + xi_sdf;
	N_sdf.rbc =  @(v) diff(v,1)+v.*xi_val(domain_x(end));
end  
     
N_sdf.op =  @(v) aa2.*diff(v,2)+bb2.*diff(v,1)+cc2.*v;
	tic;
	u0 = chebfun(@(x) 1, domain_x);
	t = [0 0.0001 0.01 0.5 1 1.5 2 2.5 3 4 5 6 7 8 9 10 15 20 25 30 50];
	prefs = cheboppref();
	prefs.discretization = @chebcolloc2;
	u = expm(N_sdf, t, u0, prefs);
	nx = 1000;
	xx = linspace(domain_x(1), domain_x(end),nx);
	for i=1:size(t,2)
		pr_el{i} = sh_el{i} -(diff(u{i},1)./u{i}.*sigmacheb+alphacheb_val);%#ok<AGROW>
		pr_el_xx(i,:) = feval(pr_el{i},xx);%#ok<AGROW>
	end
	toc;
	figure;
	plot(t,pr_el_xx(:,int16(nx/25)), 'linewidth', 2);hold on;
	plot(t,pr_el_xx(:,int16(nx/2)), 'linewidth', 2); 
	plot(t,pr_el_xx(:,int16(24*nx/25)), 'linewidth', 2); 
	ylim([min(min(pr_el_xx))-0.01 max(max(pr_el_xx))+0.01]);
	xlabel('time, years'); ylabel(m_name);title(strcat(sh_type, ', type 1 elasticity, price exposure'));
	legend(num2str(xx(int16(nx/25)), '%4.2f'),num2str(xx(int16(nx/2)), '%4.2f'),...
		num2str(xx(int16(24*nx/25)), '%4.2f'), 'Location', 'best')
	
	tic;
	u0 = chebfun(@(x) 1, domain_x);
	t = [0 0.0001 0.01 0.5 1 1.5 2 2.5 3 4 5 6 7 8 9 10 15 20 25 30 50];
	prefs = cheboppref();
	prefs.discretization = @chebcolloc2;
	u = expm(N_sdf, t, u0, prefs);
	u02 = chebfun(@(x) alphacheb_val(x), domain_x);
	u2 = expm(N_sdf, t, u02, prefs);
	nx = 1000;
	xx = linspace(domain_x(1), domain_x(end),nx);
	for i=1:size(t,2)
		pr_el2{i} = sh_el2{i}-u2{i}./u{i};%#ok<AGROW>
		pr_el_xx2(i,:) = feval(pr_el2{i},xx);%#ok<AGROW>
	end
	toc;
	figure;
	plot(t,pr_el_xx2(:,int16(nx/25)), 'linewidth', 2);hold on;
	plot(t,pr_el_xx2(:,int16(nx/2)), 'linewidth', 2); 
	plot(t,pr_el_xx2(:,int16(24*nx/25)), 'linewidth', 2); 
	ylim([min(min(pr_el_xx2))-0.01 max(max(pr_el_xx2))+0.01]);
	xlabel('time, years'); ylabel(m_name);title(strcat(sh_type, ', type 2 elasticity, price exposure'));
	legend(num2str(xx(int16(nx/25)), '%4.2f'),num2str(xx(int16(nx/2)), '%4.2f'),...
		num2str(xx(int16(24*nx/25)), '%4.2f'), 'Location', 'best')

