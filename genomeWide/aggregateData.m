function [obsBis, obsOx, C] = aggregateData(T, arg)
%aggregate CpGs. The arg argument defines the aggregated regions either 
%by a chosen window size (given in base pairs) or by containing the vector
%of split points -- HAS TO BE CORRECTED TO COPE WITH THE DIFFERENT
%CHROMOSOMES (IF WINDOW SIZE IS 1 IT DOES NOT MATTER)


%get matrix of observations
obsMatrix = table2array(T(:,2:end));
%remove all-nan rows of obsMatrix
obsMatrix = obsMatrix(~all(isnan(obsMatrix),2),:);

numOfCpGs = size(obsMatrix, 1);

firstCpGPosition = obsMatrix(1,1);
lastCpGPosition = obsMatrix(end,1);
%check the second argument
if isscalar(arg);
    windowSize = arg;
    edges = [firstCpGPosition:windowSize:lastCpGPosition, lastCpGPosition];
else
    edges = arg;
end

%MAYBE THIS CAN BE IMPLEMENTED FASTER EVENTHOUGH IT IS ALREADY FAST
%assign CpGs to genome windows and then to bins (ordered version of the filled windows)
CpGsToWindows = discretize(obsMatrix(:,1), edges);
[~, ~, CpGsToBins] = unique(CpGsToWindows);
bins = unique(CpGsToBins);
numOfBins = size(bins, 1);
obsAggr = zeros(numOfBins, size(obsMatrix, 2)-1);

%construct the cell array. 1st col is the b number
%and the second contains the CpGs' positions of the bin as string
% numOfCpGs = size(CpGsToBins, 1);
C = cell(numOfBins, 2);
C(:,1) = num2cell(bins);
%get the array of the aggregated data (1st col is the bin number)
for i=1:numOfCpGs
    %add the CpG data to
    obsAggr(CpGsToBins(i), :) = obsAggr(CpGsToBins(i), :) + nansum(obsMatrix(i, 2:end), 1);
    %put CpG position in the row of the corresponding bin (window)
    C{CpGsToBins(i), 2} = [C{CpGsToBins(i), 2}, obsMatrix(i,1)];
end

%both BS,oxBS
obsBisDay0 = obsAggr(:, 1:4);
obsBisDay3 = obsAggr(:, 9:12);
obsBisDay6 = obsAggr(:, 17:20);
obsOxDay0 = obsAggr(:, 5:8);
obsOxDay3 = obsAggr(:, 13:16);
obsOxDay6 = obsAggr(:, 21:24);

%TET_KO
% obsBisDay0 = obsAggr(:, 1:4);
% obsBisDay3 = obsAggr(:, 5:8);
% obsBisDay6 = obsAggr(:, 9:12);

obsBis = zeros(3, 4, numOfBins);
obsOx = zeros(3, 4, numOfBins);

for i=1:numOfBins
    obsBis(:,:,i) = [obsBisDay0(i,:); obsBisDay3(i,:); obsBisDay6(i,:)];
    obsOx(:,:,i) = [obsOxDay0(i,:); obsOxDay3(i,:); obsOxDay6(i,:)];
end



end