function DB_Maha = DaviesBouldinMahalanobis(numOfClusters, labels, centroids, pointMinDistances, psi)

S = zeros(1, numOfClusters);
for i = 1:numOfClusters
    S(i) = mean(sqrt(pointMinDistances(labels==i)));
end

M = zeros(numOfClusters, numOfClusters);
R = zeros(numOfClusters, numOfClusters);
D = zeros(numOfClusters,1);
for i=1:numOfClusters
    for j=1:numOfClusters
        %check here whether we should take into account the symmetric term
        %as well
        if j~=i
            M(i,j) = sqrt((centroids(:,i) - centroids(:,j))' * ((psi(:,:,i) + psi(:,:,j)) \ (centroids(:,i) - centroids(:,j))));
            R(i,j) = (S(i) + S(j)) / M(i,j);
            
        end    
    end
    D(i) = max(R(i,:));
end    

DB_Maha = mean(D);

end