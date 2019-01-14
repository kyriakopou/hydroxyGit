function [ obsBis, obsOx ] = normalizeData( obsBis, obsOx, numOfCpGs )
%NORMALIZEDATA Summary of this function goes here
%   Detailed explanation goes here
smallestSumBis = min(sum(obsBis, 2),[], 3);
smallestSumOx = min(sum(obsOx, 2), [], 3);

freqBis = obsBis ./ repmat(sum(obsBis, 2), 1, 4, 1);
freqOx = obsOx ./ repmat(sum(obsOx, 2), 1, 4, 1);

obsBis = ceil(freqBis .* repmat(smallestSumBis, 1, 4 , numOfCpGs));
obsOx = ceil(freqOx .* repmat(smallestSumOx, 1, 4 , numOfCpGs));


end

