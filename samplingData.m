function counts = samplingData(file)

numOfSamples = 500;
numObsStates = 4;

T = readtable(file, 'Delimiter',',');
%Matrices that contain all the observations
%taken from bisulfite and Oxbisulfite
obs = table2array(T(:, 3:numObsStates+2));

%get discrete prob distribution of the data for each dataPoint
pData = obs ./ repmat(sum(obs, 2), 1, numObsStates);
cumDist = cumsum(pData, 2);

counts = zeros(4, 4);
for k=1:numOfSamples
    %get 4 random numbers
    r = rand(4,1);
    for i=1:4
        counts(i, find(cumDist(i,:) >= r(i), 1)) = counts(i, find(cumDist(i,:) >= r(i), 1)) + 1;
    end
end

end