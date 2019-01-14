function [similarity, falsePosRate, truePosRate] = lucasIdea()
%LUCASIDEA Summary of this function goes here
%   Detailed explanation goes here

%read all regions
[~, ~, ~, obsBis_IAP, ~] = ... 
    UIreadDataNew('IAP_BS.txt');
[~, ~, ~, obsOx_IAP, ~] = ... 
    UIreadDataNew('IAP_oxBS.txt');
[~, ~, ~, obsBis_L1mdA, ~] = ... 
    UIreadDataNew('L1mdA_BS.txt');
[~, ~, ~, obsOx_L1mdA, ~] = ... 
    UIreadDataNew('L1mdA_oxBS.txt');
[~, ~, ~, obsBis_L1mdT, ~] = ... 
    UIreadDataNew('L1mdT_BS.txt');
[~, ~, ~, obsOx_L1mdT, ~] = ... 
    UIreadDataNew('L1mdT_oxBS.txt');
[~, ~, ~, obsBis_mSat, ~] = ... 
    UIreadDataNew('mSat_BS.txt');
[~, ~, ~, obsOx_mSat, ~] = ... 
    UIreadDataNew('mSat_oxBS.txt');
[~, ~, ~, obsBis_MuERVL, ~] = ... 
    UIreadDataNew('MuERVL_BS.txt');
[~, ~, ~, obsOx_MuERVL, ~] = ... 
    UIreadDataNew('MuERVL_oxBS.txt');
[~, ~, ~, obsBis_Afp, ~] = ... 
    UIreadDataNew('Afp_BS.txt');
[~, ~, ~, obsOx_Afp, ~] = ... 
    UIreadDataNew('Afp_oxBS.txt');
[~, ~, ~, obsBis_Zim3, ~] = ... 
    UIreadDataNew('Zim3_BS.txt');
[~, ~, ~, obsOx_Zim3, ~] = ... 
    UIreadDataNew('Zim3_oxBS.txt');
[~, ~, ~, obsBis_Ttc25, ~] = ... 
    UIreadDataNew('Ttc25_BS.txt');
[~, ~, ~, obsOx_Ttc25, ~] = ... 
    UIreadDataNew('Ttc25_oxBS.txt');
[~, ~, ~, obsBis_Snrpn, ~] = ... 
    UIreadDataNew('Snrpn_BS.txt');
[~, ~, ~, obsOx_Snrpn, ~] = ... 
    UIreadDataNew('Snrpn_oxBS.txt');


% '/Users/kyriakopou/Documents/Hydroxymethylation/dataFilesNew/newFormat/'

%concatenate all observations
obsBis = cat(3, obsBis_IAP, obsBis_L1mdA, obsBis_L1mdT, obsBis_mSat, obsBis_MuERVL, obsBis_Afp, obsBis_Zim3, obsBis_Ttc25, obsBis_Snrpn);
obsOx = cat(3, obsOx_IAP, obsOx_L1mdA, obsOx_L1mdT, obsOx_mSat, obsOx_MuERVL, obsOx_Afp, obsOx_Zim3, obsOx_Ttc25, obsOx_Snrpn);

trueSplitPoints = [5, 18, 23, 26, 30, 35, 43, 49]; 
numOfDays = 4;
numOfCpGs = size(obsBis, 3);

%normalize all regions to have approximately
%the same number of sample wrt to the smallest number of
%sample per line (OPTIONAL)
[obsBis, obsOx] = normalizeData(obsBis, obsOx, numOfCpGs);


%probability of data per day
probBisDayCpG = zeros(numOfDays, numOfCpGs);
probOxDayCpG = zeros(numOfDays, numOfCpGs);
probBisDayCpGFromGamma = zeros(numOfDays, numOfCpGs);
probOxDayCpGFromGamma = zeros(numOfDays, numOfCpGs);
obsBisRegion = zeros(numOfDays, size(obsBis, 2));
obsOxRegion = zeros(numOfDays, size(obsOx, 2));
similarity = zeros(1, numOfCpGs);
similarityPrevCpG = zeros(1, numOfCpGs);
similarityThisCpG = zeros(1, numOfCpGs);

%n is the total number of observations per day,
%X is the data matrix, A contains is the dirichlet params vector 
%for each day


%get the dirichlet params initially from the 1st CpG
A_bis = obsBis(:,:,1) + 1;
A_ox = obsOx(:,:,1) + 1; 
predSplitPoints = [];

for i=5:numOfCpGs;
     
%     if ismember(i, 49:numOfCpGs)
%        obsBis(:,:,i) = obsBis(:,:,i) .* 10;
%        obsOx(:,:,i) = obsOx(:,:,i) .* 10;
%     end    
    
    %observations of the current CpG
    X_bis = obsBis(:,:,i);
    X_ox = obsOx(:,:,i);
    
      
    %compute similarity as the logProb (1 side)
    similarityThisCpG(i) = logProbDataGivenDirichlet(X_bis, X_ox, A_bis, A_ox, numOfDays);
    if i>1
        similarityPrevCpG(i) = logProbDataGivenDirichlet(obsBis(:,:,i-1), obsOx(:,:,i-1), X_bis+1, X_ox+1, numOfDays);
    else
        similarityPrevCpG(i) = similarityThisCpG(i);
    end
    %similarity(i) = max(similarityPrevCpG(i), similarityPrevCpG(i));
    similarity(i) = (similarityThisCpG(i) + similarityPrevCpG(i)) / 2;
    
    %if the CpG does not fit to the current model update the params
    %--> we have found a new region.
    if (i>1 & similarity(i) < 2*similarity(i-1))
        predSplitPoints = [predSplitPoints, i];
        obsBisRegion = obsBis(:,:,i);
        obsOxRegion = obsOx(:,:,i);
    %else add the CpG data to the data of the region     
    else
        obsBisRegion = obsBisRegion + obsBis(:,:,i);
        obsOxRegion = obsOxRegion + obsOx(:,:,i);
        
    end
    
    %update the dirichlet params considering data of all the clustered CpGs
    A_bis = obsBisRegion + 1;
    A_ox = obsOxRegion + 1; 

    
    %check if the means of the observations are close to the data for which we estimated params 
%     meanObsBis = repmat(nBis, 1, 4) .* A_bis ./ repmat(a0_bis, 1, 4);
%     meanObsOx = repmat(nOx, 1, 4) .* A_ox ./ repmat(a0_ox, 1, 4);
    
end

%plot similarity, trueSplitPoints and predSplitPoints
%plot(1:numOfCpGs, similarityThisCpG, 1:numOfCpGs, similarityPrevCpG);
% plot(1:numOfCpGs, similarityThisCpG);
plot(1:numOfCpGs, similarity); 
hold on;
% legend('thisCpG', 'prevCpg', 'Location', 'southeast');
% legend('thisCpG', 'Location', 'southeast');
y1=get(gca,'ylim');
for x1 = trueSplitPoints
    hold on;
    plot([x1 x1], y1);
end
for x2 = predSplitPoints
    hold on;
    plot([x2 x2], y1, '--', 'LineWidth', 3);
end

%compute true positive rate (sensitivity)
%and false positive rate (1-specificity)
P = size(trueSplitPoints, 2);
N = numOfCpGs - P;
TP = size(intersect(predSplitPoints, trueSplitPoints), 2);
FP = size(predSplitPoints, 2) - TP;

truePosRate = TP / P;
falsePosRate = FP / N;



end

