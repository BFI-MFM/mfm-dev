function expectation = Expectation(den,integrand)

dx = den.x(2)-den.x(1);
expectation = sum(integrand(den.x).*den.d.*dx);

