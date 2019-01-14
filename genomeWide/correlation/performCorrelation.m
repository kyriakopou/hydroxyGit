function performCorrelation(matfileName, matCompactFileName, matClusterFileName, day, data)

%load output files 
T_full = loadTable(matfileName); 
T_compact = loadTable(matCompactFileName);

%get only rows with 3 obs days (the content of this function has to be changed if we run tet ko experiment)
[~, ~, threeDaysRows] = getRowsOfTable(T_full);

%keep only three days CpGs from T_BI
T = T_compact(threeDaysRows,:);
% T = T_compact(threeDaysRows,:);

%read T_clusters
% T_clusters = loadTable(matClusterFileName);

%clusterRows
% cluster1Rows = T_clusters.Cluster == 1;
% cluster2Rows = T_clusters.Cluster == 2;


lag = 5;
maxLag = 10000;

%create T_day for computing the correlations of a certain day
if strcmp(day, 'd0') 
    if strcmp(data, 'both')
        T_day = T(:, [1,2,5:8,17:19]);
    else
        T_day = T(:, [1,2,17:19]);
    end
elseif strcmp(day, 'd4')
    if strcmp(data, 'both')
        T_day = T(:, [1,2,9:12,20:22]);
    else
        T_day = T(:, [1,2,20:22]);
    end
else 
    if strcmp(data, 'both')    
        T_day = T(:, [1,2,13:16,23:25]);
    else
        T_day = T(:, [1,2,23:25]);
    end
end    
        
writeCorrelationToFile(T_day, 0:lag:maxLag, day)

end
