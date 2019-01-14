function [avgBisSamplesPerRegion, avgOxSamplesPerRegion, avgSamplesPerDay] = avgNumOfSamplesPerCpG()
%counts the average number of samples perDay and perCpG for the given regions
%   Detailed explanation goes here
regions = {'L1mdT', 'L1mdA', 'mSat', 'IAP', 'Afp', 'DMR6', 'DMR11', 'MuERVL'};
CpGsCell = {[1:5], [1:13], [1:3], [2 3 4 6], [1:5], [1:6], [1:8], [1:4]};

avgBisSamplesPerRegion = zeros(4, size(regions, 2));
avgOxSamplesPerRegion = zeros(4, size(regions, 2));

for i=1:size(regions,2)
    
    
    for j=CpGsCell{i}    
    
        fileName = strcat('dataFilesNew/singleCpGs/', regions{i}, 'CpG', int2str(j));
        [~, ~, ~, obsBis, obsOx, ~, ~, ~] = readData(fileName);

        avgBisSamplesPerRegion(:,i) = sum(obsBis, 2) + avgBisSamplesPerRegion(:,i); 
        avgOxSamplesPerRegion(:,i) = sum(obsOx, 2) + avgOxSamplesPerRegion(:,i);
    end
    
    avgBisSamplesPerRegion(:,i) = avgBisSamplesPerRegion(:,i) ./ size(CpGsCell{i}, 2);
    avgOxSamplesPerRegion(:,i) = avgOxSamplesPerRegion(:,i) ./ size(CpGsCell{i}, 2);
      
end

avgBisSamplesPerDay = mean(mean(avgBisSamplesPerRegion));
avgOxSamplesPerDay = mean(mean(avgOxSamplesPerRegion));
avgSamplesPerDay = (avgBisSamplesPerDay + avgOxSamplesPerDay) / 2;


