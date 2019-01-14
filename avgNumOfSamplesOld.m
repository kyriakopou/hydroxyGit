function [avgBisSamplesPerDay, avgOxSamplesPerDay, avgSamplesPerDay] = avgNumOfSamples()
%count the average number of samples perDay for the given regions
%   Detailed explanation goes here
regions = {'L1mdT', 'L1mdA', 'mSat', 'IAP', 'Afp', 'MuERVL', 'DMR6', 'DMR11'};

avgBisSamplesPerDay = zeros(4,1);
avgOxSamplesPerDay = zeros(4,1);

for i = 1:size(regions, 2);
    fileName = strcat('dataFilesNew/Wt/', regions{i});
    [~, ~, ~, obsBis, obsOx, ~, ~, ~] = readData(fileName);

    avgBisSamplesPerDay = sum(obsBis, 2) + avgBisSamplesPerDay; 
    avgOxSamplesPerDay = sum(obsOx, 2) + avgOxSamplesPerDay;

end

avgBisSamplesPerDay = mean(avgBisSamplesPerDay / size(regions, 2));
avgOxSamplesPerDay = mean(avgOxSamplesPerDay / size(regions, 2));
avgSamplesPerDay = (avgBisSamplesPerDay + avgOxSamplesPerDay) / 2;


