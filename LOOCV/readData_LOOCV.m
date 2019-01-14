function [dataPointsName, dataPoints, labels, obsBisTrain, obsOxTrain, obsBisTest, obsOxTest, errorsBisTrain, errorsOxTrain, process] = readData_LOOCV(dataPath, trainCpG, testCpG)
%read Data from the files Bis.txt and Ox.txt  

%define the number of observable states
numObsStates = 4;

%initialization
obsBisTrain = zeros(4, 4);
obsOxTrain = zeros(4, 4);
errorsBisTrain = 0;
errorsOxTrain = 0;

%read the trainCpGs and merge them in one table
for i = trainCpG
    %read the observations and the error constants from
    %the appropriate .txt file and store in tables
    file1 = strcat(dataPath, 'CpG', int2str(i), '_BS.txt');
    file2 = strcat(dataPath, 'CpG', int2str(i), '_oxBS.txt');

    Tbis = readtable(file1, 'Delimiter', ',');
    Tox = readtable(file2, 'Delimiter', ',');   

    %Matrices that contain all the observations
    %taken from bisulfite and Oxbisulfite
    obsBis = table2array(Tbis(:, 3:numObsStates+2));
    obsOx = table2array(Tox(:, 3:numObsStates+2));
    %sum them up
    obsBisTrain = obsBisTrain + obsBis;
    obsOxTrain = obsOxTrain + obsOx;
    
    %column vectors of the conversion errors for
    %each time point
    errorsBis = table2array(Tbis(:, numObsStates+3:numObsStates+5));
    errorsOx = table2array(Tox(:, numObsStates+3:numObsStates+5));
    errorsBisTrain = errorsBisTrain + errorsBis;
    errorsOxTrain = errorsOxTrain + errorsOx;
    
    
    %vector containing the days of the observations
    dataPoints = table2array(Tbis(:, 1));
    %array containing the labels of the dataPoints
    labels = table2array(Tbis(:, 2));
    %last column defines if the transition btw two dataPoints 
    %is a cell-replication or not
    process = table2array(Tbis(:, numObsStates+6));
    %first variable name (either day or dataPoint)
    dataPointsName = Tbis.Properties.VariableNames(1);
    
end

errorsBisTrain = errorsBisTrain ./ size(trainCpG, 2);
errorsOxTrain = errorsOxTrain ./ size(trainCpG, 2);


%read the testCpG data
file1 = strcat(dataPath, 'CpG', int2str(testCpG), '_BS.txt');
file2 = strcat(dataPath, 'CpG', int2str(testCpG), '_oxBS.txt');

Tbis = readtable(file1, 'Delimiter', ',');
Tox = readtable(file2, 'Delimiter', ',');   

%Matrices that contain all the observations
%taken from bisulfite and Oxbisulfite
obsBisTest = table2array(Tbis(:, 3:numObsStates+2));
obsOxTest = table2array(Tox(:, 3:numObsStates+2));



end