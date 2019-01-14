%store results in the output table (all chromosomes) -- vertical
%concatenation of each of the previous tables
T_resultsMLE_WG = vertcat(T_resultsMLEChr{:});
T_resultsBI_WG = vertcat(T_resultsBIChr{:});
MLEFilePath = strcat(resultsPath, '/T_resultsMLE_WG');
BIFilePath = strcat(resultsPath, '/T_resultsBI_WG');
save(MLEFilePath, 'T_resultsMLE_WG');
save(BIFilePath, 'T_resultsBI_WG');

%get more compact versions of T_resultsMLE_WG and T_resultsBI_WG
% getCompactTable(T_resultsMLE_WG);
% getCompactTable(T_resultsBI_WG);

 %write results of MLE table to a .txt file
txt_MLEFilePath = createTxtFile(T_resultsMLE_WG, MLEFilePath);
txt_BIFilePath = createTxtFile(T_resultsBI_WG, BIFilePath);

igvMLEFile = copyTxtToIgv(txt_MLEFilePath);
igvBIFile = copyTxtToIgv(txt_BIFilePath);