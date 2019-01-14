function numResultsRow = getDetailedNumerResults(params, Cov, pAllStates, KLBis, KLOx, hyperVol, accRate)

switch nargin 
    case 6
        accRate = nan;
end   
%get the upper triangle's elements of Cov matrix
CovElements = Cov(triu(true(size(Cov))));

%store results into table's row
numResultsRow = [reshape(pAllStates', [1 27]), params, CovElements', KLBis, KLOx, hyperVol, accRate];
        

end