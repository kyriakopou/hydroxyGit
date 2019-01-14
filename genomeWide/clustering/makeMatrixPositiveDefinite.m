function A = makeMatrixPositiveDefinite(A)

A = topdm(A);
while ~isPositiveDefinite(A)
    A = nearestSPD(A);
end

end