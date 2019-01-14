function SSE = elbow(numOfClusters, labels, pointMinDistances)

S = zeros(1, numOfClusters);
for i = 1:numOfClusters
    S(i) = sum(pointMinDistances(labels==i));
end
SSE = sum(S);

end