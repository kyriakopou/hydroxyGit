function DB = DaviesBouldinEuclid(numOfClusters, labels, centroids, pointMinDistances)

S = zeros(1, numOfClusters);
for i = 1:numOfClusters
%     S(i) = sqrt(mean(sum(pointMinDistances(labels==i))));
    %average euclidean distance of all points in cluster i from each
    %centroid
    S(i) = mean(sqrt(pointMinDistances(labels==i)));
end

M = zeros(numOfClusters, numOfClusters);
R = zeros(numOfClusters, numOfClusters);
D = zeros(numOfClusters,1);
for i=1:numOfClusters
    for j=1:numOfClusters
        %check here whether we should also take into account the symmetric
        %term
        if j~=i
            M(i,j) = sqrt((centroids(:,i) - centroids(:,j))' * (centroids(:,i) - centroids(:,j))); 
            R(i,j) = (S(i) + S(j)) / M(i,j);
        end    
    end
    D(i) = max(R(i,:));
end    

DB = mean(D);

end