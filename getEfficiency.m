%returns the efficiency values for time points of vector t
function y = getEfficiency(coef, t)

    t = t';
    numOfTimePoints = size(t, 1);
    %GETEFFICIENCY Summary of this function goes here
    %   Detailed explanation goes here
    numOfCoef = size(coef, 2);
    tPowers = zeros(numOfTimePoints, numOfCoef);
    for i=1:numOfCoef
        tPowers(:,i) = t'.^(i-1);
    end

    y = sum(repmat(coef, [numOfTimePoints,1]) .* tPowers, 2); 

end

