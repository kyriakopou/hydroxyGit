function b = isPositiveDefinite(M)

b = all(eig(M) > 0);

end