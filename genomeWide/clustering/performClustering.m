function T_clusters = performClustering(matfileName, matCompactFileName, clusteringMethod, ...
    dataDesc, distMeasure, criterion)
%CLUSTERINGEFFICIENCIES Summary of this function goes here
%   Detailed explanation goes here
%The function reads from the WT or TET results files and clusters every CpG
%wrt the enzymes efficiencies, methylation levels, or both

if nargin == 5
    % criterion = 'CalinskiHarabasz';
    criterion = 'DaviesBouldin';
end

% distMeasure = 'cosine';

%load output file 
T_final = loadTable(matfileName); 
T_compact = loadTable(matCompactFileName);

%get only rows with 3 obs days
[~, ~, threeDaysRows] = getRowsOfTable(T_final);

%keep only three days CpGs from T_BI
T_final = T_final(threeDaysRows,:);
T_compact = T_compact(threeDaysRows,:);

%compress further (TEST CASE)
%get random tables rows (CpGs)
% randRows = randi(size(T_final, 1), [10^5, 1]);
% randRows = sort(randRows);
% T_final = T_final(randRows,:);
% T_compact = T_compact(randRows,:);

% T_input = T_input(threeDaysRows,:);
%get the number of observations for each CpG
% numOfObsPerCpG = sum(T_input{:,3:end}, 2);

eff = zeros(size(T_final, 1), 2);



%---NEW .MAT FORMAT 
% %maintenance
eff(:,1) = T_final.maint_b0;
eff(:,2) = T_final.maint_b1;

%deNovo
eff(:,3) = T_final.deNovo_b0;
eff(:,4) = T_final.deNovo_b1;

%hydroxy
eff(:,5) = T_final.hydroxy_b0;
eff(:,6) = T_final.hydroxy_b1;

eff(:,7) = T_final.pRecogn;


if strcmp(dataDesc, 'both')
    X = [eff, levels];
elseif strcmp(dataDesc, 'eff') 
    X = eff;
else
    X = levels;
end


CovEl = T_final{:,39:66};
numOfParams = 7;
numOfElements = size(X,1);


%get covariances matrices from the upper diagonal elements
Covariances = nan(numOfParams, numOfParams, numOfElements);
if strcmp(clusteringMethod, 'kerror')
    for i=1:numOfElements
        Covariances(:,:,i) = getCovarianceFromUpperDiag(CovEl(i,:));    
    end
%     
%     Covariances = eye(numOfParams, numOfParams);
%     Covariances = repmat(Covariances, [1 1 numOfElements]);

    Covariances = massageCovs(Covariances);
    %make invOfCovs positive-semi definite
    invOfCovs = multinv(Covariances);

    %make invOfCovs symmetric 
    invOfCovs = (invOfCovs + permute(invOfCovs, [2 1 3])) / 2;
    invOfCovs = massageCovs(invOfCovs);
end


%%choose the type of the clustering
maxNumOfClusters = 10;
numOfReps = 10;
numOfClusters = 1:maxNumOfClusters;
DB = zeros(1, length(numOfClusters));
CH = zeros(1, length(numOfClusters));
SSE = zeros(1, length(numOfClusters));

%do k-means euclidean (matlab)
% eva = clustering(X(:,1:6), clusteringMethod, numOfClusters, distMeasure, criterion);
% plot(eva);

DBmin = inf;
totAvg = mean(X(:,1:6), 1)';
for m = numOfClusters
    if strcmp(clusteringMethod, 'kmeans')
        [idx, centroids, minDist] = kmeansMod(X(:,1:6)', m, numOfReps);        
        DB(numOfClusters == m) = DaviesBouldinEuclid(m, idx, centroids, minDist);
        CH(numOfClusters == m) = CalinskiHarabaszEuclid(m, idx, centroids, totAvg, minDist);
        SSE(numOfClusters == m) = elbow(m, idx, minDist);
        
    else        
        [idx, centroids, minDist, psi, totAvg, psiTotAvg] = kError(X(:,1:6)', m, invOfCovs(1:6,1:6,:), numOfReps);
        DB(numOfClusters == m) = DaviesBouldinMahalanobis(m, idx, centroids, minDist, psi);
        CH(numOfClusters == m) = CalinskiHarabaszMahalanobis(m, idx, centroids, totAvg, minDist, psi, psiTotAvg);
        SSE(numOfClusters == m) = elbow(m, idx, minDist);
    end     
    
    %keep the solution with the smallest DB so far
    if m > 1
        if DB(numOfClusters == m) < DBmin 
            DBmin = DB(numOfClusters == m);
            optNumOfClusters = m;
            idxOpt = idx;
            centroidsOpt = centroids;
            minDistOpt = minDist;
        end
    end
end


%optimal number of clusters
fprintf('The optimal number of clusters is %d \n', optNumOfClusters);
idxOpt = idxOpt';


%create table T_cluster and write it to a .txt file
T_clusters = createClusteringTable(idxOpt, T_final);
outputPath = strcat('~/Desktop/DeepResults/WT/clustering/', clusteringMethod, '/', distMeasure, '_', ...
    num2str(optNumOfClusters), 'clusters');

pdfFileName_db = strcat(outputPath, '_db.pdf');
pdfFileName_ch = strcat(outputPath, '_ch.pdf');
pdfFileName_elbow = strcat(outputPath, '_elbow.pdf');
pdfFileName_clusters = strcat(outputPath, '_clusters.pdf');

%save sse and db arrays
txtFileName = strcat(outputPath, '.txt');
dlmwrite(strcat(outputPath, '_db.txt'), DB);
dlmwrite(strcat(outputPath, '_ch.txt'), CH);
dlmwrite(strcat(outputPath, '_elbow.txt'), SSE);

matOutFileName = strcat(outputPath, '.mat');
save(matOutFileName, 'T_clusters');
writetable(T_clusters, txtFileName, 'Delimiter', '\t');

%change the file's extension from txt to .igv
% copyTxtToIgv(txtFileName);

fig1 = figure();
plot(numOfClusters(2:end), DB(2:end));
xlabel('Number of clusters');
ylabel('Davies-Bouldin value');
saveas(fig1, pdfFileName_db);

fig2 = figure();
plot(numOfClusters, CH);
xlabel('Number of clusters');
ylabel('Calinski Harabasz');
saveas(fig2, pdfFileName_ch);

fig3 = figure();
plot(numOfClusters, SSE);
xlabel('Number of clusters');
ylabel('Within clusters sum of squares');
saveas(fig3, pdfFileName_elbow);

%plot the clusters
fig4 = figure();
plotClusters(T_clusters, T_compact, dataDesc);
saveas(fig4, pdfFileName_clusters);

end


