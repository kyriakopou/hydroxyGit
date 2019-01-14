function Cov = makeCovSymmetric(Cov)

if (~issymmetric(Cov))
    Cov = triu(Cov) + triu(Cov, +1)';
end

end