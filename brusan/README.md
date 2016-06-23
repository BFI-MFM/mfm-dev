“A Macroeconomic Model with a Financial Sector”  by
Markus K. Brunnermeier and Yuliy Sannikov

The main program is solve_equilibrium: to compute equilibrium type solve_equilibrium.m in the
command window and press enter. The file will then produce two figures:

-Figure 1 shows important equilibrium quantities such as q, psi, drift and volatility of eta.
-Figure 2 shows expert and household utility within the model.

To compute the equilibrium under a different set of parameters change line 14 of
solve_equilibrium.m, in which parameter values are assigned. Investment adjustment cost
parameter can be modified directly in the function investment.m.

- fnct.m is a key function that programs the set of equations to compute the derivatives
of and using the equations from Proposition 2 of the paper.

- solve_equilibrium.m operates by searching for initial conditions near   eta=0

- evntfcn.m is used to determine when integration of the differential equations should be
terminated for a given set of initial conditions, and how the initial conditions should be
modified on the next iteration.