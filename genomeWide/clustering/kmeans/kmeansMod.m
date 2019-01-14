function [labelMin, muMin, sqrdEuclidMin] = kmeansMod(X, m, numOfReps)
% Perform kmeans clustering.
% Input:
%   X: d x n data matrix
%   m: initialization parameter
% Output:
%   label: 1 x n sample labels
%   mu: d x k center of clusters
%   energy: optimization target value
% Written by Mo Chen (sth4nth@gmail.com).

%precompute X'X
dotX = dot(X,X,1);
objFuncMin = inf;

for rep=1:numOfReps

    label = init(X, m);
    n = numel(label);
    idx = 1:n;
    last = zeros(1,n);
    
    %the algorithm takes advantage of the fact that min(X'X + mu'mu -
    %2X'mu) is the same as distance = min(mu'mu/2 - mu'X). To get the euclidean
%     distance of each point u need sqrdEuclid = X'X + 2*distance
    iter = 0;
    while any(label ~= last)
        iter = iter + 1;
        [~,~,last(:)] = unique(label);                  % remove empty clusters
        mu = X*normalize(sparse(idx,last,1),1);         % compute cluster centers 
        [distance, label] = min(dot(mu,mu,1)'/2-mu'*X,[],1);  % assign sample labels
%         [euclidDist, label] = min(pdist2(X', mu'), [], 2);
    end
    
    sqrdEuclid = dotX + 2*distance;
    objFunc = sum(sqrdEuclid);
    
    %print info for this repetition
    fprintf('rep %d: iterations = %d : ObjFunction = %.4f \n\n', rep, iter, objFunc);
    
    if objFunc < objFuncMin
        objFuncMin = objFunc;
        sqrdEuclidMin = sqrdEuclid;
        labelMin = label;
        muMin = mu;
%         euclidDistMin = euclidDist; 
    end

    
end

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

