function nSPDRows = findNonSpdRows(A)

n = size(A, 3);
nSPDRows = false(1, n);
for i=1:n
    if ~isPositiveDefinite(A(:,:,i))
        nSPDRows(i) = true;
    end
end


end