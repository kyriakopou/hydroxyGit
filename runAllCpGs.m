function runAllCpGs()

% regions = {'L1mdT', 'L1mdA', 'mSat', 'IAP', 'Afp', 'DMR6', 'DMR11', 'MuERVL'};
% CpGsCell = {[1:5], [1:13], [1:3], [2 3 4 6], [1:5], [1:6], [1:8], [1:4]};

% regions = {'Snrpn'};
% CpGsCell = {[1:6]};

regions = {'Ttc25', 'Zim3'};
CpGsCell = {[1:6], [1:8]};

for i=1:size(regions,2)
    for j=CpGsCell{i}
        fileName = strcat('dataFilesNew/singleCpGs/', regions{i}, 'CpG', int2str(j));
        [x, sigma] = estimateDSHydroxy(fileName, 'oxBS');
        
    end
end

end