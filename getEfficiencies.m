function [EFF, varEFF] = getEfficiencies(x, maxDegree, timePoints, varargin)
%given the parameter vector containing the coefficients of rates as
%polynomial bases functions of maxDegree compute at each timepoint (of vector timePoints) 
%the efficiencies together with their variances

switch nargin
    case 3
        knots = [];
        Cov = [];
    case 4
        knots = varargin{1};
        Cov = [];
    case 5
        knots = varargin{1};
        Cov = varargin{2};
end

numOfKnots = length(knots);
numOfCoefProRate = maxDegree+1+numOfKnots;
numOfTimePoints = length(timePoints);
numOfRates = 5;

%powers of the polynomial plus one term for each knot
powers = [0:maxDegree, repmat(maxDegree, [1 numOfKnots])];
powersMat = repmat(powers', [1, numOfTimePoints]);

%timePoints for the polynomial terms 
timePointsMat = repmat(timePoints, [numOfCoefProRate-numOfKnots, 1]);
%check if knots is empty
if ~isempty(knots) 
    knotsMat = repmat(knots', [1, numOfTimePoints]);
    timePointsWithKnotsMat = repmat(timePoints, [numOfKnots, 1]) > knotsMat  .* (repmat(timePoints, [numOfKnots, 1]) - knotsMat);
else
    timePointsWithKnotsMat = [];
end    

%one term timePoint-knot for each knot
timePointsMat = [timePointsMat; timePointsWithKnotsMat];

%each row tPowers contains all terms of the quadratic spline for a specific time point
tPowers = (timePointsMat .^ powersMat)';

paramsMatMu1 = repmat(x(1:numOfCoefProRate), [numOfTimePoints, 1]);
mu1 = sum(paramsMatMu1 .* tPowers, 2);

paramsMatMu2 = repmat(x(numOfCoefProRate+1:2*numOfCoefProRate), [numOfTimePoints, 1]);
mu2 = sum(paramsMatMu2 .* tPowers, 2);

paramsMatMuh = repmat(x(2*numOfCoefProRate+1:3*numOfCoefProRate), [numOfTimePoints, 1]);
h = sum(paramsMatMuh .* tPowers, 2);

paramsMatMuf = repmat(x(3*numOfCoefProRate+1:4*numOfCoefProRate), [numOfTimePoints, 1]);
f = sum(paramsMatMuf .* tPowers, 2);

paramsMatMudem = repmat(x(4*numOfCoefProRate+1:5*numOfCoefProRate), [numOfTimePoints, 1]);
dem = sum(paramsMatMudem .* tPowers, 2);

EFF = [mu1, mu2, h, f, dem];

%return also the variances of the efficiencies over time for all
%timepoint and all efficiencies
if nargin == 5
    varEFF = zeros(numOfTimePoints, numOfRates);

    for i = 1:numOfRates    
     
        rateInd = (i-1)*numOfCoefProRate;

        coefCovMatrix = zeros(numOfCoefProRate, numOfCoefProRate, numOfTimePoints);
        for k=1:numOfTimePoints
            CovMatrixOfRate = Cov(rateInd+1:rateInd+numOfCoefProRate, rateInd+1:rateInd+numOfCoefProRate);   
            coefCovMatrix(:,:,k) = tPowers(k,:)' * tPowers(k,:); 
        end
        
        varEFF(:,i) = sum(sum(coefCovMatrix .* repmat(CovMatrixOfRate, [1, 1, numOfTimePoints]))); 
    
    end
    
end



end