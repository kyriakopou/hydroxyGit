function [model, params, maxDegree, knots, matrixFileName] = produceMatricesFiles(fileName)

% hidVarValues = getVarValues;

%READ THE MODEL FILES AND PRODUCE THE DISCRETE TIME TRANSITION MATRICES
[model, hiddenVars, params, maxDegree, knots, processes, rulesLHS, rulesRHS, probs] = ... 
    txtInputFileParser(fileName);

numOfParams = length(params);
numOfProcesses = size(processes, 2);
%number of CpG sites
numOfSites = getNumOfSites;
%create StateSet
hidVarValues = hiddenVars;
hidStateSet = createStateSet(hidVarValues, numOfSites);
numOfHidStates = size(hidStateSet, 2);

numOfCoef = maxDegree + 1 + size(knots, 2);
numOfRates = (numOfParams -1) / numOfCoef;

matrixFileName = cell(1, numOfProcesses+1);


%for each process that we read from .txt file produce
for proc=1:numOfProcesses
    
    %Check if fileName exists. Otherwise produce the file
    a = strsplit(num2str(knots));
    matrixFileName{proc} = strcat(processes{proc}, '_deg_', int2str(maxDegree));
    if ~isempty(knots)
        matrixFileName{proc} = strcat(matrixFileName{proc}, '_', a{:});
    end    
    fileName = strcat('~/Documents/code/MATLAB/Hydroxymethylation/modelGeneration/matrices/process/active/', ...
        matrixFileName{proc}, '.m');
    
    if ~exist(fileName, 'file')
        
        fprintf('Producing necessary discrete time matrixFile for %s. It might take a while... \n', processes{proc});    

        %symbolic matrices
        matrixSym = sym(diag(ones(1, numOfHidStates)));
        derSym = sym(zeros(numOfHidStates, numOfHidStates, numOfParams));
        secDerSym = sym(zeros(numOfHidStates, numOfHidStates, numOfParams, numOfParams));


        %for each state get all descendents and fill the three symbolic matrices
        %(transMatrix, firstDer, secDer)
        for i=1:numOfHidStates
            [descendants, transProbs] = rulesToState(getStateFromIndex(i), rulesLHS{proc}, rulesRHS{proc}, probs{proc});
            for j=1:size(descendants, 2)
                %index of the descendant
                indJ = getIndexFromState(descendants{j});
                %fill the entry i,indJ of each of the matrices
                matrixSym(i,indJ) = transProbs(j);

                %find the variables of transition probability
                vars = symvar(transProbs(j));
                derProb = sym(zeros(1, size(vars, 2)));
                secDerProb = sym(zeros(1, size(vars, 2)));

                for k=1:size(vars, 2)
                    %get the derivative of transition prob
                    derProb(k) = diff(transProbs(j), vars(k));
                    secDerProb(k) = diff(derProb(k), vars(k));
                    %find the place of vars(k) in params vector
                    pInd = find(params == vars(k));
                    %place it in the correct derMatrix
                    derSym(i, indJ, pInd) = derProb(k);
                    secDerSym(i, indJ, pInd, pInd) = secDerProb(k);

                end    

            end
        end

        %fill the rest derivative and second derivative matrices
        %MAYBE MORE EFFICIENT S.T. NOT SO MANY SYMBOLIC OPERATIONS
        syms t;

        polyTerms = t.^(0:maxDegree);
        truncTerms = heaviside(t-knots) .* (t - knots).^(maxDegree);
        baseTerms = [polyTerms, truncTerms];
        baseTermsSecDer = reshape(baseTerms, [numOfCoef 1]) * baseTerms;

        for i=1:numOfRates

            rateInd = (i-1)*numOfCoef;
            %---FIRSTDERIVATIVE---
            for j1=1:numOfCoef
                derSym(:,:,rateInd+j1) = derSym(:,:,rateInd+1) * baseTerms(j1);
                %----SECDERIVATIVE-----
                for j2=1:numOfCoef
                    secDerSym(:,:,rateInd+j1,rateInd+j2) = baseTermsSecDer(j1,j2) * secDerSym(:,:,rateInd+1,rateInd+1);
                end
            end
        end
        
        matlabFunction(matrixSym, derSym, secDerSym, 'File', fileName, 'Outputs', {'M', 'derM', 'secDerM'});
        disp('done!')
    end
    
end





end