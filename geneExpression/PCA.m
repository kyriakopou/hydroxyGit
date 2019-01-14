function [varExplained, avgCorrTrain, avgCorrTest] = PCA(fileName, numOfRepeats)
%Get a table of features and compute the 
%number of components (linear comb) needed to acchieve
%some given variance level
[X, y, headers] = readAndreasFile(fileName);

[PC, SCORES, latent] = pca(X);

varExplained = cumsum(latent)./sum(latent);
numOfPCToGet95Var = find(varExplained >= 0.95, 1);

maxNumOfPCs = size(X,2);
plot(1:maxNumOfPCs, 100*varExplained, '-bo');
xlabel('numOfComponents');
ylabel('Percent of Variance explained');

% SCORES(:,1:2)
figure();
biplot(PC(:,1:2), 'varlabels', headers(1:end-1));

%%------REGRESSION-----------
%%
totNumOfLines = size(SCORES, 1);
lines = (1:totNumOfLines)';
numOfTrainLines = floor(70/100 * totNumOfLines);

avgCorrTrain = zeros(1, maxNumOfPCs);
avgCorrTest = zeros(1, maxNumOfPCs);
stdCorrTrain = zeros(1, maxNumOfPCs);
stdCorrTest = zeros(1, maxNumOfPCs);

for numOfComps=1:maxNumOfPCs    
    
    corrTrain = zeros(1, numOfRepeats);
    corrTest = zeros(1, numOfRepeats);
    
    for k=1:numOfRepeats
        
        trainLines = randperm(totNumOfLines, numOfTrainLines)';
        testLines = setdiff(lines, trainLines);

        Z_train = SCORES(trainLines,1:numOfComps);
        Z_test = SCORES(testLines, 1:numOfComps);
        y_train = y(trainLines);
        y_test = y(testLines);

        %train data
        b = regress(y_train, Z_train);
        
        %see pred for train data
        pred_y = Z_train*b;
%         scatter(pred_y, y_train);
        corrCoefMatTest = corrcoef(pred_y, y_train);
        corrTrain(k) = corrCoefMatTest(1,2);

        %see pred for test data    
        pred_y = Z_test*b;
%         scatter(pred_y, y_test);
        corrCoefMatTest = corrcoef(pred_y, y_test);
        corrTest(k) = corrCoefMatTest(1,2);
        
    end    
    
    avgCorrTrain(numOfComps) = mean(corrTrain);
    stdCorrTrain(numOfComps) = std(corrTrain);
    avgCorrTest(numOfComps) = mean(corrTest);
    stdCorrTest(numOfComps) = std(corrTest);
    
end


figure();
errorbar(1:maxNumOfPCs, avgCorrTrain, 1.96*stdCorrTrain);
legend('')
hold on;
errorbar(1:maxNumOfPCs, avgCorrTest, 1.96*stdCorrTest);
legend('trainCorr','testCorr');

xlabel('numOfComponents')
ylabel('corr (pred vs real gene expr.)');
    


end


