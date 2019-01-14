function CH = CalinskiHarabaszEuclid(numOfClusters, labels, centroids, totMean, pointMinDistances)

%within cluster variance for each cluster
SSw = zeros(1, numOfClusters);
SSb = zeros(1, numOfClusters);
for i=1:numOfClusters
    SSw(i) = sum(pointMinDistances(labels == i));
    n_i = length(pointMinDistances(labels == i));
    SSb(i)  = n_i * ((centroids(:,i) - totMean)' * (centroids(:,i) - totMean));
end

%overall within cluster variance
SSw = sum(SSw);
%between cluster variance
SSb = sum(SSb);

N = length(labels);
CH = SSb/SSw * (N - numOfClusters)/(numOfClusters-1);

end