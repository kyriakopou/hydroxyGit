function [avgCorrCoefTrain, avgCorrCoefTest] = PLSRegression(fileName, numOfRepeats)

[X, y, headers] = readAndreasFile(fileName);
totNumOfLines = size(X,1);
numOfPredictors = size(X,2);
lines = (1:totNumOfLines)';
numOfTrainLines = floor(70/100 * totNumOfLines);
maxNumOfComponents = numOfPredictors;


%train and test corr matrices
avgCorrCoefTrain = zeros(1, maxNumOfComponents);
stdCorrCoefTrain = zeros(1, maxNumOfComponents);
avgCorrCoefTest = zeros(1, maxNumOfComponents);
stdCorrCoefTest = zeros(1, maxNumOfComponents);

for numOfComponents = 1:maxNumOfComponents
    
    corrCoefTest = zeros(1, numOfRepeats);
    corrCoefTrain = zeros(1, numOfRepeats);
    for k = 1:numOfRepeats
        
        trainLines = randperm(totNumOfLines, numOfTrainLines)';
        testLines = setdiff(lines, trainLines);

        % numOfTestLines = size(X,1) - numOfTrainLines;
        X_train = X(trainLines,:);
        X_test = X(testLines, :);
        y_train = y(trainLines);
        y_test = y(testLines);

        %do PLS regression
        [~, theta_train, Z_train, ~, ~, ~, ~, stats] = plsregress(X_train, y_train, numOfComponents);

        % X0 = X - repmat(mean(X), [size(X, 1), 1]);
        y0 = y - mean(y);
        y0_train = y0(trainLines);
        y0_test = y0(testLines);
        phi_train = stats.W;

        %train correlation error
%         figure();
%         scatter(Z*theta', y0_train);
        corrTrain = corrcoef(Z_train * theta_train', y0_train);
        corrCoefTrain(k) = corrTrain(1,2);

        %test correlation error
        Z_test = X_test * phi_train;
        y_predTest = Z_test * theta_train';
%         figure();
%         scatter(y_predTest, y0_test);
        corrTest = corrcoef(y_predTest, y0_test);
        corrCoefTest(k) = corrTest(1,2);
        
    end    
    
    avgCorrCoefTrain(numOfComponents) = mean(corrCoefTrain);
    stdCorrCoefTrain(numOfComponents) = std(corrCoefTrain);
        
    avgCorrCoefTest(numOfComponents) = mean(corrCoefTest);
    stdCorrCoefTest(numOfComponents) = std(corrCoefTest);  
    
end   
figure();
errorbar(1:maxNumOfComponents, avgCorrCoefTrain, 1.96*stdCorrCoefTrain);
legend('')
hold on;
errorbar(1:maxNumOfComponents, avgCorrCoefTest, 1.96*stdCorrCoefTest);
legend('trainCorr','testCorr');

xlabel('numOfComponents')
ylabel('corr (pred vs real gene expr.)');
    
end