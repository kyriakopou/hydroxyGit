function [xCorrDelta, P, RL, RU] = getEffCorr(T, delta)

%get pos+delta and pos-delta vectors to find all CpG pairs of distance delta
varNames = {'Chromosome', 'Start'};

T_subPlusDelta = table(T.Chromosome, T.Start + delta, 'VariableNames', varNames);
T_subMinusDelta = table(T.Chromosome, T.Start - delta, 'VariableNames', varNames);

[~, ~, idxIntoPosFromPlus] = intersect(T_subPlusDelta(:,1:2), T(:,1:2));
[~, ~, idxIntoPosFromMinus] = intersect(T_subMinusDelta(:,1:2), T(:,1:2));


% n_delta = length(idxIntoPosFromPlus);
% muPlus = mean(T{idxIntoPosFromPlus, 3:end});
% muMinus = mean(T{idxIntoPosFromMinus, 3:end});
% stdPlus = std(T{idxIntoPosFromPlus, 3:end});
% stdMinus = std(T{idxIntoPosFromMinus, 3:end});
% 
% auto_delta = 1./((n_delta-1).*stdPlus.*stdMinus) .* sum((T{idxIntoPosFromPlus, 3:end} - muPlus) .*  ...
%     (T{idxIntoPosFromMinus, 3:end} - muMinus));

%get all the columns apart from the positions
A = [T{idxIntoPosFromPlus, 3:end}, T{idxIntoPosFromMinus, 3:end}];
[xCorrDelta, P, RL, RU] = corrcoef(A);
%get only the correlation between minus and plus CpGs
s = size(xCorrDelta, 1);
xCorrDelta = xCorrDelta(1:s/2, s/2+1:end); 
P = P(1:s/2, s/2+1:end);
RL = RL(1:s/2, s/2+1:end);
RU = RU(1:s/2, s/2+1:end);

end