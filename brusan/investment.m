function [Phi iota] = investment(q)
% [Phi iota] = investment(q)
% finds the optimal investment rate iota and rate of capital creation Phi
% given the price q

theta = 10;  % adjustment cost parameter
Phi = (q - 1)/theta;  iota = Phi + theta*Phi^2/2;
