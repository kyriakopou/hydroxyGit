function [pData, pModel, pAllStates ] = UIestimateLevels(dataPoints, obsMatrix, errorsMatrix)


numOfDays = size(dataPoints, 1);


%number of observable states
numObsStates = size(obsMatrix, 2);
numOfExperiments = sum(any(any(obsMatrix)));
%numOfStates depending on what experim. we did
numOfHiddenStates = getNumberOfStatesTemp(numOfExperiments);

%throw out experiments with all zero observations from obsMatrix
indToRem = find(sum(any(obsMatrix), 2) == 0);
obsMatrix(:,:,indToRem) = [];
errorsMatrix(:,:,indToRem) = [];



EMatrix = zeros(numOfHiddenStates, numObsStates, numOfDays, numOfExperiments);
% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
for t=1:numOfDays
    EMatrix(:,:,t,1) = myConvErrorsBisAct(1-errorsMatrix(t,:,1));
    if numOfExperiments > 1 
        EMatrix(:,:,t,2) = myConvErrorsOxAct(1-errorsMatrix(t,:,2));
    end    
    if numOfExperiments > 2    
        EMatrix(:,:,t,3) = myConvErrorsMSIAct(1-errorsMatrix(t,:,3));
    end    
end


%using MLEstimator
pData = zeros(numOfDays, numObsStates, numOfExperiments);
pModel = zeros(numOfDays, numObsStates, numOfExperiments);
pAllStates = zeros(numOfDays, numOfHiddenStates);

for t=1:numOfDays
    E = EMatrix(:,:,t,:);
    EMatrixForTimeT = reshape(E, [size(E, 1), size(E,2), size(E, 4), size(E, 3)]);
    [pAllStates(t,:), pData(t,:,:), pModel(t,:,:)] = maxLikelihood(obsMatrix(t,:,:), EMatrixForTimeT); 
end


