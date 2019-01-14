function [avgCorrTrain, avgCorrTest] = svmRegression(fileName)

[X, y, ~] = readAndreasFile(fileName);

totNumOfLines = size(X, 1);
maxNumOfVars = size(X, 2);
lines = (1:totNumOfLines)';
numOfTrainLines = floor(70/100 * totNumOfLines);

avgCorrTrain = zeros(1, maxNumOfVars);
avgCorrTest = zeros(1, maxNumOfVars);


MdLin = fitrsvm(X, y, 'KFold', 10);
% MdPoly = fitrsvm(X, y, 'KFold', 10, 'KernelFunction', 'polynomial');
% MdGaussian = fitrsvm(X, y, 'KFold', 10, 'KernelFunction', 'gaussian');




end

