function igvFile = copyTxtToIgv(txtFile)

%change the BI file's extension to .igv
igvFile = strcat(txtFile(1:end-4), '.igv');
copyfile(txtFile, igvFile);


end