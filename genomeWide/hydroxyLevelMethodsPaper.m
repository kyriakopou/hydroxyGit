function hydroxyLevelMethodsPaper(regionName)
%GIVEPASCAL Summary of this function goes here
%   Detailed explanation goes here

fileName_BS = strcat(regionName, '_BS.txt');
fileName_oxBS = strcat(regionName, '_oxBS.txt');
fileName_errors = strcat(regionName, '_errors.txt');
[~, ~, obsBis, ~, ~] = UIreadDataNew(fileName_BS);
[~, ~, obsOx, ~, ~] = UIreadDataNew(fileName_oxBS);
[~, dataPoints, dataLabels, errorsBis, errorsOx, ~] = UIreadErrors(fileName_errors);


numOfDays = size(dataPoints, 1);
numObsStates = size(obsBis, 2);
numOfStates = 9;
dataType = 'oxBS';

% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
EBis = zeros(numOfStates, numObsStates, numOfDays);
EOx = zeros(numOfStates, numObsStates, numOfDays);
for t=1:numOfDays
    EBis(:,:,t) = convErrorBis(errorsBis(t,:));
    EOx(:,:,t) = convErrorOx(errorsOx(t,:));   
end


pHidFromData = zeros(numOfDays, numOfStates);
pBis = zeros(numOfDays, numObsStates);
pOx = zeros(numOfDays, numObsStates);
pBisModel = zeros(numOfDays, numObsStates);
pOxModel = zeros(numOfDays, numObsStates);

for t=1:numOfDays;
    [pHidFromData(t,:), pBis(t,:), pOx(t,:), pBisModel(t,:), pOxModel(t,:)] = maxLikelihood(obsBis(t,:), obsOx(t,:), EBis(:,:,t), EOx(:,:,t), dataType); 

    uu = pHidFromData(t,1);
    um = sum(pHidFromData(t,2:3), 2);
    uh = sum(pHidFromData(t,4:5), 2);
    hm = sum(pHidFromData(t,6:7), 2);
    hh = pHidFromData(t,9);
    toth = hm + uh + hh;
    mm = pHidFromData(t,8);
    hydroxyStates = [uh, hm, hh];
    plotedStates = [mm, toth, um, uu];
    %createFigurePaper(plotedStates, hydroxyStates, days, regionName);
    pAllStates = pHidFromData;

end            

%THIS PATH HAS TO BE GIVEN AS AN ARGUMENT BY THE USER!
outFilePath = strcat('/Users/kyriakopou/Desktop/pascal/wholeGenome/results/', regionName, '.txt');

fileID = fopen(outFilePath, 'w');

fprintf('\n');
disp('Optimization is done! Check the file  ')
disp(outFilePath);
fprintf('\n');


fprintf(fileID, '---------------------------------------------------------------\n');

C = num2cell(pBis);
C = [dataLabels, C];
formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

fprintf(fileID, 'The data distribution of BS states \n\n');
fprintf(fileID, 'data TT TC CT CC \n');
[nrows, ncols] = size(C);
for row = 1:nrows
    fprintf(fileID, formatSpec,C{row,:});
end
fprintf(fileID, '\n');

C = num2cell(pBisModel);
C = [dataLabels, C];
formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

fprintf(fileID, 'The predicted distribution of BS states \n\n');
fprintf(fileID, 'data TT TC CT CC \n');
[nrows, ncols] = size(C);
for row = 1:nrows
    fprintf(fileID, formatSpec,C{row,:});
end

fprintf(fileID, '---------------------------------------------------------------\n');

C = num2cell(pOx);
C = [dataLabels, C];
formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

fprintf(fileID, 'The data distribution of oxBS states \n\n');
fprintf(fileID, 'data TT TC CT CC \n');
[nrows, ncols] = size(C);
for row = 1:nrows
    fprintf(fileID, formatSpec,C{row,:});
end
fprintf(fileID, '\n');

C = num2cell(pOxModel);
C = [dataLabels, C];
formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

fprintf(fileID, 'The predicted distribution of oxBS states \n\n');
fprintf(fileID, 'data TT TC CT CC \n');
[nrows, ncols] = size(C);
for row = 1:nrows
    fprintf(fileID, formatSpec,C{row,:});
end

fprintf(fileID, '---------------------------------------------------------------\n');

C = num2cell(pAllStates);
C = [dataLabels, C];
formatSpec = '%s %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n';

fprintf(fileID, 'The predicted distribution of the hidden states \n\n');
fprintf(fileID, 'data uu um mu uh hu hm mh mm hh \n');
[nrows, ncols] = size(C);
for row = 1:nrows
    fprintf(fileID, formatSpec,C{row,:});
end
fprintf(fileID, '\n');
    


end

