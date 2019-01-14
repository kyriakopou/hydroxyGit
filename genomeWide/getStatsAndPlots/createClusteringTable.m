function T_clusters = createClusteringTable(idx, T_final)
%create the clustering table and copy its contents to a .txt with txtFileName and to an igv file

%create empty feature column
Feature = zeros(size(T_final, 1), 1);
%store cluster for every CpG in T_out and copy this to .mat, .txt, .igv file 
varNames = {'Chromosome', 'Start', 'End', 'Methylation', 'Cluster'};
T_clusters = table(T_final.Chromosome, T_final.Start, T_final.End, Feature, idx, 'VariableNames', varNames);
% save(strcat(outputFilePath, 'clustering/T_clusteringBI_', num2str(optNumOfClusters), '_', criterion), 'T_clusters');



end