function eva = clustering(X, method, numOfClusters, distMeasure, criterion)
%CLUSTERINGEFFICIENCIES Summary of this function goes here

% myFunc = @(X,K)(clusterdata(X, K, 'distance', 'seuclidean'));

clusteringFun = str2func(method);

if strcmp(method, 'kmeans')
    % evaluate the clusterings for a specific range of clusters
    myFunc = @(X,K)(clusteringFun(X, K, 'distance', distMeasure, 'emptyaction', 'drop',...
         'replicates', 1, 'Display', 'final'));
end

% evaluate what is the best number of clusters according to C.H. criterion
eva = evalclusters(X, myFunc, criterion, 'KList', numOfClusters);


end

