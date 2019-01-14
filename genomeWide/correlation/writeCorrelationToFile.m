function writeCorrelationToFile(T_day, lags, day)

numOfRVs = size(T_day(:,3:end), 2) ;
xCorr = zeros(numOfRVs, numOfRVs, length(lags));
P = zeros(numOfRVs, numOfRVs, length(lags));
RL = zeros(numOfRVs, numOfRVs, length(lags));
RU = zeros(numOfRVs, numOfRVs, length(lags));

for i=lags
    [xCorr(:,:,lags==i), P(:,:,lags==i), RL(:,:,lags==i), RU(:,:,lags==i)] = getEffCorr(T_day, i);      
end


%store the lags and the correlations in a file
outputPath = strcat('~/Desktop/DeepResults/WT/correlation/');
corrTxtFileName = strcat(outputPath, 'corr_', day, '.txt');
pStatTxtFileName = strcat(outputPath, 'pStat_', day, '.txt');
lowerBoundTxtFileName = strcat(outputPath, 'LB_', day, '.txt');
upperBoundTxtFileName = strcat(outputPath, 'UB_', day, '.txt');

dlmwrite(corrTxtFileName, xCorr);
dlmwrite(pStatTxtFileName, P);
dlmwrite(lowerBoundTxtFileName, RL);
dlmwrite(upperBoundTxtFileName, RU);




end