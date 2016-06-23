function [value,isterminal,direction] = evntfcn(eta, F) 
% [value,isterminal,direction] = evntfcn(eta, F) returns three
% variables used by ode45 to determine when to terminate integration
%
% written by Yuliy Sannikov

global qmax
value = [(qmax - F(3)) F(2) F(4)];  % difference between qmax and q, 
                                    % first derivative of theta, first derivative of q
isterminal = [1 1 1];   % terminate computation in all three cases
direction = [0 0 0];    % event occurs whether we get there from above or below 
                                    
