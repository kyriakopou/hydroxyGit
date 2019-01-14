function [avgCorrTrain, avgCorrTest] = linRegression(fileName, numOfRepeats)

[X, y, ~] = readAndreasFile(fileName);

totNumOfLines = size(X, 1);
maxNumOfVars = size(X, 2);
lines = (1:totNumOfLines)';
numOfTrainLines = floor(70/100 * totNumOfLines);

avgCorrTrain = zeros(1, maxNumOfVars);
avgCorrTest = zeros(1, maxNumOfVars);
stdCorrTrain = zeros(1, maxNumOfVars);
stdCorrTest = zeros(1, maxNumOfVars);

for numOfVars=1:maxNumOfVars    
    
    corrTrain = zeros(1, numOfRepeats);
    corrTest = zeros(1, numOfRepeats);
    
    for k=1:numOfRepeats
        
        trainLines = randperm(totNumOfLines, numOfTrainLines)';
%         trainLines = 1:numOfTrainLines';
        testLines = setdiff(lines, trainLines);

        X_train = X(trainLines, 1:numOfVars);
        X_test = X(testLines, 1:numOfVars);
        y_train = y(trainLines);
        y_test = y(testLines);

        %train data
        b = regress(y_train, X_train);
        
        %see pred for train data
        pred_y = X_train*b;
%         scatter(pred_y, y_train);
        corrCoefMatTest = corrcoef(pred_y, y_train);
        corrTrain(k) = corrCoefMatTest(1,2);

        %see pred for test data    
        pred_y = X_test*b;
%       scatter(pred_y, y_test);
        corrCoefMatTest = corrcoef(pred_y, y_test);
        corrTest(k) = corrCoefMatTest(1,2);
        
    end    
    
    avgCorrTrain(numOfVars) = mean(corrTrain);
    stdCorrTrain(numOfVars) = std(corrTrain);
    avgCorrTest(numOfVars) = mean(corrTest);
    stdCorrTest(numOfVars) = std(corrTest);
    
end


figure();
errorbar(1:maxNumOfVars, avgCorrTrain, 1.96*stdCorrTrain);
legend('')
hold on;
errorbar(1:maxNumOfVars, avgCorrTest, 1.96*stdCorrTest);
legend('trainCorr','testCorr');

xlabel('numOfVars')
ylabel('corr (pred vs real gene expr.)');

end