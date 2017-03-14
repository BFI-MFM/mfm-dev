function [ergodic_density,  p_X] = ergodic_density_mc(model_name)
%computes density for a macro model using Monte-Carlo simulations
%input: model name
%output: ergodic density and bins
load (model_name);

d = [domain_x(1), domain_x(2)];

T =10000;
T_stat = 9900;
dt = 0.1;
dist_percentile =[0.2500 0.5000 0.7500 ];
N=20000;

X = zeros(N,T);
X(:,1) = linspace(d(1),d(2),N);

for t = 2 : T
    mu_t = mucheb(X(:,t-1));
    sigma_t = sigmacheb(X(:,t-1));
    
    Xtemp = X(:,t-1) + mu_t * dt + sigma_t .* randn(N,1) * dt^0.5;
    
    % treatment of boundaries
    % left boundary
    if strcmp(l.boundary_type,'natural')
        %Neumann - natural
        Xtemp(Xtemp<d(1)) = d(1);
        X(:,t) = Xtemp;
    elseif strcmp(l.boundary_type, 'reflecting')
        Xtemp(Xtemp<d(1)) = d(1);
       % Xtemp(Xtemp>d(1)) = 2*d(1) - Xtemp(Xtemp>d(1));
        X(:,t) = Xtemp;
    end
    
    % right boundary
    if strcmp(r.boundary_type,'natural')
        %Neumann - natural
        Xtemp(Xtemp>d(2)) = d(2);
        X(:,t) = Xtemp;
    elseif strcmp(r.boundary_type, 'reflecting')
         Xtemp(Xtemp>d(2)) = d(2);
        %Xtemp(Xtemp>d(2)) = 2*d(2) - Xtemp(Xtemp>d(2));
        X(:,t) = Xtemp;
    end
end

X = X(:,1+T_stat:end);
x_size = size(X(:),1);

X = sort(X(:));
state_percentile = X(min(x_size,max(1,round(dist_percentile*x_size))))';

[p_v, p_X] = hist(X,200);
ergodic_density = p_v/x_size/(p_X(2)-p_X(1));

zz=max(ergodic_density);
lines_p= linspace(0,zz,100);
figure;plot(p_X, ergodic_density, 'k', ...
    state_percentile(1).*ones(1,100),lines_p,'--r', ...
    state_percentile(2).*ones(1,100),lines_p,'--b',...
    state_percentile(3).*ones(1,100),lines_p,'--g',...
    'Linewidth',1.6)
legend('density','25%', '50%','75%', 'Location', 'best');