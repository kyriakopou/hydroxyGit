function [labelMin, muMin, minDistMin, psi, muMahaData, psiMahaData] = kError(X, numOfClusters, invOfCovs, numOfReps)
% Perform kErrors clustering as described in Kummar et al.
% Input:
%   X: d x n data matrix
%   m: initialization parameter
% Output:
%   label: 1 x n sample labels
%   mu: d x k center of clusters
%   energy: optimization target value
% Written by Mo Chen (sth4nth@gmail.com).


dim = size(X,1);


%CHECK IF SIGMA SHOULD CHANGE BEC OF THE DIM OF X
sumOfInv = zeros([size(invOfCovs(:,:,1)), numOfClusters]);
sumOfInvTimesX = zeros(dim, numOfClusters); 
mu = zeros(dim, numOfClusters);

objFuncMin = inf; 
for rep = 1:numOfReps
    lastObjFunc = inf;
    
    %random initialization of labels
%     label = init(X, numOfClusters);
    %kmeans++ initialization
    mu = initCentroids(dim, numOfClusters);
    label = init(X, mu);
      
    n = numel(label);
    last = zeros(1,n);
    distance = inf(n, numOfClusters);
    objFunc = inf;
    [minDist, ~] = min(distance, [], 2);
    
    iter = 0;
    while any(label ~= last) && iter < 101
        iter = iter + 1; 
        [~,~,last(:)] = unique(label);                  % remove empty clusters
        OldObjFunc = objFunc;
        OldMinDist = minDist;
                          
        %compute the mahalanobis cluster centroids
        for cluster=1:numOfClusters
            indOfm = (label == cluster);
            sumOfInv(:,:,cluster) = sum(invOfCovs(:,:,indOfm), 3);
            sumOfInvTimesX(:,cluster) =  sum(mtimesx(invOfCovs(:,:,indOfm), permute(X(:,indOfm), [1 3 2])), 3);
            if ~isPositiveDefinite(sumOfInv(:,:,cluster))
                sumOfInv(:,:,cluster) = makeMatrixPositiveDefinite(sumOfInv(:,:,cluster));
            end    
            
            %compute the new cluster centroid
            mu(:, cluster) = sumOfInv(:,:,cluster) \ sumOfInvTimesX(:,cluster) ;  %inv(sumOfInv(:,:,m)) * sumOfInvTimesX

            %compute the new distance of each point to the new cluster centroid
            firstProd = mtimesx(permute((X - repmat(mu(:,cluster), 1, n))', [3 2 1]), invOfCovs);
            distance(:,cluster) = mtimesx(firstProd, permute((X - repmat(mu(:,cluster), 1, n)), [1 3 2]));
        end
        
        %get new clustering (labels and minDistances)
        [minDist, label] = min(distance, [], 2);
        
        %compute the objective function value
        objFunc = sum(minDist);
        if objFunc > OldObjFunc && ~isequal(label',last)
            ind = find(label == last');
            minDist(ind) = OldMinDist(ind);
            fprintf('%f, %f \n', OldObjFunc, objFunc);
        end
        
        label = label';
        
        fprintf('iter %d: objFunc = %f \n', iter, objFunc)
    end
    
    %print info for this repetition
    fprintf('rep %d: iterations = %d : ObjFunction = %.4f \n\n', rep, iter, objFunc);

    %set the min objective function and the covariance matrices 
    %of the mahalanobis centroid of each cluster 
    if objFunc < objFuncMin
        objFuncMin = objFunc;
        minDistMin = minDist; 
        labelMin = label;
        muMin = mu;
        psi = multinv(sumOfInv);
    end    
    
    
end

%compute the mahaDataCentroid
sumOfInv = sum(invOfCovs, 3);
sumOfInvTimesX = sum(mtimesx(invOfCovs, permute(X, [1 3 2])), 3);
if ~isPositiveDefinite(sumOfInv)
    sumOfInv = makeMatrixPositiveDefinite(sumOfInv);
end    

%compute the new cluster centroid and its psi
muMahaData = sumOfInv \ sumOfInvTimesX;  %inv(sumOfInv(:,:,m)) * sumOfInvTimesX
psiMahaData = inv(sumOfInv);


%random initialization function
function label = init(X, m)
[d,n] = size(X);
if numel(m) == 1                           % random initialization
    mu = X(:,randperm(n,m));
    [~,label] = min(dot(mu,mu,1)'/2-mu'*X,[],1); 
elseif all(size(m) == [1,n])               % init with labels
    label = m;
elseif size(m,1) == d                      % init with seeds (centers)
    [~,label] = min(dot(m,m,1)'/2-m'*X,[],1); 
end

