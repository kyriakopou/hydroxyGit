function [pdf, pdfVector] = myTruncNormalPdf(x, mean, sigmas, lastDay)

numOfParams = size(x, 2);
truncRegion = zeros(1,numOfParams);
pdfVector = normpdf(x, mean, sigmas);

constPos = 1:2:numOfParams;
gradPos = 2:2:numOfParams;
numOfConstParams = size(constPos, 2);

%get probability mass of the truncated regions for both const and gradient
%params
% truncRegion(constPos) = normcdf(ones(1, numOfConstParams), mean(constPos), sigmas(constPos)) - ...
%     normcdf(zeros(1, numOfConstParams), mean(constPos), sigmas(constPos));
% truncRegion(gradPos) = normcdf((1-x(constPos(1:end-1)))./lastDay, mean(gradPos), sigmas(gradPos)) -  ...
%     normcdf(-x(constPos(1:end-1))./lastDay, mean(gradPos), sigmas(gradPos));

%standardize the variables and %get probability mass of the truncated regions for both const and gradient params
truncRegion(constPos) = phi((ones(1, numOfConstParams) - mean(constPos)) ./ sigmas(constPos)) -  ...
    phi((zeros(1, numOfConstParams) - mean(constPos)) ./ sigmas(constPos));
truncRegion(gradPos) = phi(((1-x(constPos(1:end-1)))./lastDay - mean(gradPos)) ./ sigmas(gradPos)) -  ... 
    phi((-x(constPos(1:end-1))./lastDay - mean(gradPos)) ./ sigmas(gradPos));

%normalize all normals by their truncated region
pdfVector = pdfVector ./ truncRegion;
pdf = prod(pdfVector);

%warning message
% if isnan(pdf) || isinf(pdf)
%     disp('stop');
% end    
   


