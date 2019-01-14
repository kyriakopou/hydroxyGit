function logProb = dirichletMultinomialPmf(x, a)
%DIRICHLETMULTINOMIALPMF Summary of this function goes here
%   Detailed explanation goes here

n = sum(x);
% a0 = sum(a);

%compute the prob using the beta function
% prob = n*beta(n, a0) / prod(nonzeros(x * beta(x, a0)));

%compute prob using gamma function
% prob = factorial(n) * gamma(a0) / gamma(n+a0) * prod(gamma(x+a) ./ x.* gamma(a));

%logProb using gamma function
firstTwoSums = nansum(log(1:n)) - nansum(log(n+4:2*n+3));
lastTwoSums = 0;
for k = 1:4
    lastTwoSums = lastTwoSums + nansum(log(nonzeros(x(k):2*x(k)))) - nansum(log(1:x(k)));
end    

logProb = firstTwoSums + lastTwoSums;


end

