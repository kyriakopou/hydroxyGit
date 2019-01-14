function compareInputOutput(inputFilePath, outputFilePath, errorsFile)
%MODELFIT Summary of this function goes here
%   Detailed explanation goes here

T_in = loadTable(inputFilePath);
T_out = loadTable(outputFilePath);

%get value of general vars -- THIS TO BE READ (GENERALIZED)
numOfExperiments = 2;
numObsStates = 4;
numOfHiddenStates = 9;
numOfDays = 3;

%read errors and get EMatrix for the conv errors for each day and
%experiment
EMatrix = zeros(numOfHiddenStates, numObsStates, numOfDays, numOfExperiments);
[~, ~, daysLabels, errorsMatrix, ~] = UIreadErrors(errorsFile);

% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
for t=1:numOfDays
    EMatrix(:,:,t,1) = convErrorBis(errorsMatrix(t,:,1));
    if numOfExperiments > 1 
        EMatrix(:,:,t,2) = convErrorOx(errorsMatrix(t,:,2));
    end    
end


%%
%FOR THE WHOLE FIT SECTION WE HAVE TO CHANGE THE MODEL RESULTS BY
%CONSIDERING THE ERRORS (WE NEED FOR THIS ALL HIDDENS STATES - NEW RUN)
%--FIT DATA VS MODEL PLOTS-- 

[TT_bsAll, hemi_bsAll, CC_bsAll] = plotModelFitGW_KO(T_in, T_out, EMatrix, daysLabels);




end

