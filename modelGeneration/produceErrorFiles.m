function produceErrorFiles()

%initializations
numOfSites = getNumOfSites;
hidVarValues = getVarValues;
obsVarValues = {'T', 'C'};
hidStateSet = createStateSet(hidVarValues, numOfSites);
obsStateSet = createStateSet(obsVarValues, numOfSites);
numOfHidStates = size(hidStateSet, 2);
numOfObsStates = size(obsStateSet, 2);


%READ THE ERROR FILE AND PRODUCE THE ERROR MATRICES
[~, ~, ~, experiments, rulesLHS, rulesRHS, probs] = txtInputFileParser('/errorFiles/errorActive.txt');
numOfExperiments = size(experiments, 2);

%for all experiments that we read from .txt file
for exper=3:numOfExperiments

    %symbolic matrices
    matrixSym = sym(diag(ones(1, numOfObsStates)));

    %for each state get all descendents and fill the three symbolic matrices
    %(transMatrix, firstDer, secDer)
    for i=1:numOfHidStates
        state = getStateFromIndex(i);
        [descendants, transProbs] = rulesToState(state, rulesLHS{exper}, rulesRHS{exper}, probs{exper});
        for j=1:size(descendants, 2)
            %index of the descendant
            indJ = getIndexFromState(descendants{j}, obsVarValues);
            %fill the entry i,indJ of each of the matrices
            matrixSym(i,indJ) = transProbs(j);
            
        end
    end
    
    matlabFunction(matrixSym, 'File', strcat('~/Documents/code/MATLAB/Hydroxymethylation/modelGeneration/matrices/errors/', experiments{exper}), 'Outputs', {'E'});

end



end