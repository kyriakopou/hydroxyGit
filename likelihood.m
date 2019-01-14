function [logLik, derlogLik, pModel] = likelihood(p, obsMatrix, errorMatrix)

numOfStates = size(p, 2);
numObsStates = size(obsMatrix(:,:,1), 2);
numOfExperiments = size(obsMatrix, 3); 

pModel = zeros(1, numObsStates, numOfExperiments);
for i=1:numOfExperiments
    pModel(:,:,i) = p*errorMatrix(:,:,i);
end    


%Bisulfite, oxBS and msiBS normalized Observations for all time points 
n = zeros(size(obsMatrix));
for i=1:numOfExperiments
    n(:,:,i) = obsMatrix(:,:,i) ./ repmat(nansum(obsMatrix(:,:,i), 2), 1, numObsStates);
end

% nData = reshape(n, 1, i*4);

%compute likelihood summing over all states 
logLik = nansum(nansum(obsMatrix .* log(pModel)));
logLik = -logLik;
%logLik = norm(nData - pEst);

derlogLik = zeros(1, numOfStates);
for i=1:size(p, 2)
    derlogLik(i) = nansum(nansum(obsMatrix .* errorMatrix(i,:,:) ./ pModel));
end


derlogLik = -derlogLik;




end