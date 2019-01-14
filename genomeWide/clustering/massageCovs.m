function A = massageCovs(A)

nSPDRows = findNonSpdRows(A);
nRows = find(nSPDRows == 1);

for i=nRows
    A(:,:,i) = makeMatrixPositiveDefinite(A(:,:,i)); 
end


end