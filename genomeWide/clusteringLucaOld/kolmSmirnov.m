function [rejectVectorBis, rejectVectorOx] = kolmSmirnov(X1_bis, X1_ox, X2_bis, X2_ox)

prob_X1Bis = X1_bis ./ repmat(sum(X1_bis, 2), 1, 4, 1);
prob_X1Ox = X1_ox ./ repmat(sum(X1_ox, 2), 1, 4, 1);
prob_X2Bis = X2_bis ./ repmat(sum(X2_bis, 2), 1, 4, 1);
prob_X2Ox = X2_ox ./ repmat(sum(X2_ox, 2), 1, 4, 1);

distBis = max(abs(X1_bis - X2_bis), [], 2);
distOx = max(abs(X1_ox - X2_ox), [], 2);

numOfObs_X1Bis = sum(X1_bis, 2);
numOfObs_X1Ox = sum(X1_ox, 2);
numOfObs_X2Bis = sum(X2_bis, 2);
numOfObs_X2Ox = sum(X2_ox, 2);

%c for confidence level a = 0.05
c = 1.36;

rejectVectorBis = (distBis >= c * sqrt((numOfObs_X1Bis + numOfObs_X2Bis) ./ numOfObs_X1Bis .* numOfObs_X2Bis)); 
rejectVectorOx = (distOx >= c * sqrt(numOfObs_X1Ox + numOfObs_X2Ox) ./ numOfObs_X1Ox .* numOfObs_X2Ox); 



end