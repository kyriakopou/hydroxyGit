function txt_FilePath = createTxtFile(T_results_WG, outfilePath)

% %write results of BI table to a .txt file
txt_FilePath = strcat(outfilePath, 'GW_BI.txt');
writetable(T_results_WG, txt_FilePath, 'Delimiter', '\t');


end