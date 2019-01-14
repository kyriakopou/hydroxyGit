function [dataPointsName, dataPoints, labels, obsBis, obsOx, errorsBis, errorsOx, process] = readData(dataPath)
%read Data from the files Bis.txt and Ox.txt  

%define the number of observable states
numObsStates = 4;

%read the observations and the error constants from
%the appropriate .txt file and store in tables
file1 = strcat(dataPath, '_BS.txt');
file2 = strcat(dataPath, '_oxBS.txt');

Tbis = readtable(file1, 'Delimiter',',');
Tox = readtable(file2, 'Delimiter',',');

%vector containing the days of the observations
dataPoints = table2array(Tbis(:, 1));
%array containing the labels of the dataPoints
labels = table2array(Tbis(:, 2));

%Matrices that contain all the observations
%taken from bisulfite and Oxbisulfite
obsBis = table2array(Tbis(:, 3:numObsStates+2));
obsOx = table2array(Tox(:, 3:numObsStates+2));
%column vectors of the conversion errors for
%each time point
errorsBis = table2array(Tbis(:, numObsStates+3:numObsStates+5));
errorsOx = table2array(Tox(:, numObsStates+3:numObsStates+5));
%last column defines if the transition btw two dataPoints 
%is a cell-replication or not
process = table2array(Tbis(:, numObsStates+6));

%first variable name (either day or dataPoint)
dataPointsName = Tbis.Properties.VariableNames(1);


end