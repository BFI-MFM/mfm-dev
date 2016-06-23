function [fp dyn] = fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)
% [fp dyn] = fnct(eta, f, r, rho, a, a_, delta, delta_, sigma)
% takes scalar eta, a 4x1 vector f = [theta, theta', q, q'], and parameters
% r, rho, a, a_, delta, delta_, sigma, and computes fp, the derivative of f with
% respect to eta, as well as 
% dyn = [psi, sigma_eta*eta, sigma_q, mu_eta*eta, mu_q, iota, leverage, rk, r_k]
%
% assumes that the production set of experts is (a - iota, Phi(iota) - delta) and
% households, (a_ - iota, Phi(iota) - delta_), where the function [Phi iota] = investment(q)
% solves for Phi and iota
%
% written by Yuliy Sannikov

% search for psi between eta (lower bound) and min(f(3)/f(4) + eta, 1) (upper bound) 

[Phi iota] = investment(f(3));

psi_L = eta; psi_R = min(f(3)/f(4) + eta, 1);  
for n = 1:50,
            psi = (psi_L + psi_R)/2;
            amplification = 1 - f(4)/f(3)*(psi - eta);
            
            % VOLATILITY COMPUTATION
            sigma_eta_eta = sigma*(psi - eta)/amplification;  % sigma_eta *times* eta
            sigma_q = sigma_eta_eta*f(4)/f(3); 
            sigma_theta = sigma_eta_eta*f(2)/f(1);
            risk_premium = - sigma_theta*(sigma + sigma_q);
            
            household_premium = (a_ - a)/f(3) + delta - delta_ + risk_premium;
                   
            if household_premium > 0, % households want to hold more
                psi_R = psi;
            else
                psi_L = psi;
            end
end

mu_q = r - (a - iota)/f(3) - Phi + delta - sigma*sigma_q + risk_premium;
mu_eta_eta = -(psi - eta)*(sigma + sigma_q)*(sigma + sigma_q + sigma_theta) + eta*(a - iota)/f(3) + eta*(1 - psi)*(delta_ - delta);
qpp = 2*(mu_q*f(3) - f(4)*mu_eta_eta)/sigma_eta_eta^2; 
thetapp = 2*((rho - r)*f(1) - f(2)*mu_eta_eta)/sigma_eta_eta^2; 

fp = [f(2); thetapp; f(4); qpp];

leverage = psi/eta;
rk = r + risk_premium;
r_k = r + household_premium;

dyn = [psi, sigma_eta_eta, sigma_q, mu_eta_eta, mu_q, iota, leverage, rk, r_k];





