function Cov = getCovarianceFromUpperDiag(CovUpTriangle)
%get the covariance matrix from the upper diagonal elements
    numOfParams = 5;
    
    Cov = zeros(numOfParams, numOfParams);
    Cov(triu(true(numOfParams))) = CovUpTriangle;
    Cov = Cov + tril(Cov', -1);   


end