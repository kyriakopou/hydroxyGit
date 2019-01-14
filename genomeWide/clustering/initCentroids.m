function mu = initCentroids(dim, numOfClusters)

muConst = rand(dim/2, numOfClusters);
muGrad = rand(3,numOfClusters) .* (muConst/6 + muConst/6) - (muConst / 6);

mu = [muConst; muGrad];
mu = mu([1 4 2 5 3 6],:);

end