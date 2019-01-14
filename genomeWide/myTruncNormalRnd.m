function z = myTruncNormalRnd(mean, sigmas, lastDay)
%generates vector z from univariate independent truncated normal distr.
numOfParams = size(mean, 2);
%x is the vector of normalized truncated RVs
x = zeros(1, numOfParams);

lb = zeros(1, numOfParams);
ub = zeros(1, numOfParams);
%z is the non-normalized truncated vector
z = zeros(1,numOfParams);
for i=1:numOfParams
    if mod(i,2) == 1
        x(i) = trandn((0-mean(i))/sigmas(i), (1-mean(i))/sigmas(i));
        z(i) = mean(i) + sigmas(i) * x(i);
    else
        lb(i) = -z(i-1)/lastDay;
        ub(i) = (1-z(i-1))/lastDay;
        x(i) = trandn((lb(i) - mean(i))/sigmas(i), (ub(i) - mean(i)) / sigmas(i));
        z(i) = mean(i) + sigmas(i) * x(i);
    end    
end
    
end