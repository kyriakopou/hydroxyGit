function [logLik, derlogLik, pBisModel, pOxModel] = likelihoodGW(p, obsBis, obsOx, EBis, EOx)

E = [EBis, EOx]; 
pEst = p*E;
pBisModel = pEst(1:4);
pOxModel = pEst(5:8);
numOfStates = size(p, 2);
numObsStates = size(obsBis, 2);

%Bisulfite normalized Observations for all time points 
nBis = obsBis ./ repmat(nansum(obsBis, 2), 1, numObsStates);
%Oxidative Bisulfite normalized Observations for all time points
nOx = obsOx ./ repmat(nansum(obsOx, 2), 1, numObsStates);
nData = [nBis, nOx];

%compute likelihood summing over all states 
logLik = nansum(nansum(obsBis .* log(pBisModel) + obsOx .* log(pOxModel)));
logLik = -logLik;
%logLik = norm(nData - pEst);

derlogLik = zeros(1, numOfStates);
for i=1:size(p, 2)
    derlogLik(i) = nansum(obsBis .* EBis(i,:) ./ pBisModel + obsOx .* EOx(i,:) ./ pOxModel);
end
derlogLik = -derlogLik;




end