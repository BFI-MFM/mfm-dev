LW = 'linewidth'; MS = 'markersize';
disp('intersections of curves as a bivariate rootfinding problem');
disp(' here are the intersections between the splat curve and a figure-of-eight curve');
disp('Press any key to continue...') ;
pause;
tic;
t = chebfun('t',[0,2*pi]);
sp = exp(1i*t) + (1+1i)*sin(6*t).^2;     % splat curve
figof8 = cos(t) + 1i*sin(2*t);           % figure of eight curve
plot(sp,LW,1.6), hold on
plot(figof8,'r',LW,1.6), axis equal

d = [0 2*pi 0 2*pi];
f = chebfun2(@(s,t) sp(t)-figof8(s),d);  % rootfinding
r = roots(real(f),imag(f));              % calculate intersections
spr = sp(r(:,2));
plot(real(spr),imag(spr),'.k',MS,20), ylim([-1.1 2.1])
hold off
toc;

disp('FOC for 2D function');
disp('the scalar 2D function f(x,y) is represented by a chebfun2 object');
disp('a gradient of a scalar function f(x,y) is a chebfun function');
disp('here are the FOC points of a sum of  20 Gaussians:');
disp('Press any key to continue...') ;
pause;
tic;

num_gaussians = 20;
f = chebfun2(0);
lambda0 = 20;
rng('default');
for k = 1:num_gaussians;
    x0 = 2*rand-1; y0 = 2*rand-1;
    lambda = lambda0*rand+lambda0;
    f = f + chebfun2(@(x,y) exp(-lambda*((x-x0).^2 + (y-y0).^2)));
end
plot(f), hold on
r = roots(gradient(f));
plot3(r(:,1),r(:,2),f(r(:,1),r(:,2)),'k.','markersize',20)
zlim([0 4]), hold off, colormap(winter);
toc;

disp('BVP systems with unknown parameters');
disp('the scalar 2D function f(x,y) is represented by a chebfun2 object');
disp('a gradient of a scalar function f(x,y) is a chebfun function');
disp('here are the FOC points of a sum of  20 Gaussians:');
disp('Press any key to continue...') ;
pause;
tic;
omega = 1;
N = chebop(@(x, u , omega) diff(u,2) - u - sin(omega.*x/pi), [-pi pi]);
N.init = [chebfun(@(x) sin(omega*x/pi), [-pi,pi]);2];
N.lbc = @(u,omega) u - 1;
N.rbc = @(u,omega) [u - 1; diff(u) - 1];
[u1,omega] = N\0;
disp('omega');disp(omega);
plot(u1, LW, 1.6); hold on;
N.init = [chebfun(1, [-pi pi]);4];
[u2,omega] = N\0;
disp('omega');disp(omega);
plot(u2, LW, 1.6); 
legend('initial guess=sin( \omega x /\pi)','initial guess = 1','Location', 'best');hold off;
toc;