function [params, paramNames, Cov] = cleanNanParams(params, paramNames, Cov)

params(del) = [];
paramNames(del) = [];
Cov(:,del) = [];
Cov(del,:) = [];

end