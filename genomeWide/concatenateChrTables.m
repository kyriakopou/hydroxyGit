function concatenateChrTables(dirPath)

%load all chromosome files from dir
mat = dir(fullfile(dirPath, '*.mat'));
%store all chromosomes to T_BI_WG and T_MLE_WG
T_allChr = [];
for q = 1:length(mat)
    load(strcat(dirPath, mat(q).name));
    %this is needed in case we have this 60x1 cell array for each
    %chromosome
    T_resultsMLEChr(cellfun('isempty', T_resultsMLEChr)) = [];
    T_results = T_resultsMLEChr{1};
%     save(mat(q).name, 'T_results');
    
    T_allChr = vertcat(T_allChr, T_results);

end

sortrows(T_allChr, 'Chromosome', 'ascend');

%store results in the output table (all chromosomes) -- vertical
%concatenation of each of the previous tables
filePath = strcat(dirPath, 'T_MLE_WG');
save(filePath, 'T_allChr');
 
%get more compact versions of T_resultsMLE_WG and T_resultsBI_WG
% getCompactTable(T_resultsMLE_WG);
% getCompactTable(T_resultsBI_WG);


end