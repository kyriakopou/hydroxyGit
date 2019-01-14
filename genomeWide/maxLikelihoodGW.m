function [p, pbis, pox, pBisModel, pOxModel] = maxLikelihoodGW(obsBis, obsOx, EBis, EOx, dataType)
%OXBISOLVE Solve equation system to get initial distribution by estimating
%the maximum likelihood estimator for the initial probability distribution

numOfStates = 9;
numObsStates = size(obsBis, 2);

%Bisulfite normalized Observations for all time points 
pbis = obsBis ./ repmat(nansum(obsBis, 2), 1, numObsStates);
%Oxidative Bisulfite normalized Observations for all time points
pox = obsOx ./ repmat(nansum(obsOx, 2), 1, numObsStates);


% Initial guess for probability distribution
p0 = [0.5, 0.5, 0, 0, 0, 0, 0, 0, 0];
Aeq = [ones(1,9); zeros(8,9)];
beq = [1; zeros(8,1)];

%define the parameter bounds
if (strcmp(dataType, 'oxBS'))
    lb = zeros(1, numOfStates);
    ub = ones(1, numOfStates);
%bounds for david's states
% (no hydroxylated states)    
else 
    lb = zeros(1, numOfStates);
    ub = [1, 1, 1, zeros(1,4), 1, 0];
end

%FOR PARALLEL EXECUTION
% currentPool = gcp;
% if(currentPool.Connected == false)
%     parpool(4);
% end

options = optimset('Algorithm', 'interior-point', 'GradObj', 'on', 'Display', 'off', 'DerivativeCheck', 'off', 'MaxFunEvals', 100);
obj_func = @(p)likelihoodGW(p, obsBis, obsOx, EBis, EOx);
problem = createOptimProblem('fmincon', 'x0', p0, 'objective', obj_func, 'Aeq', Aeq, 'beq', beq, 'lb', lb, 'ub', ub, 'options', options);
ms = MultiStart('StartPointsToRun', 'bounds-ineqs', 'Display', 'off', 'MaxTime', 20);
p = run(ms, problem, 5);


%get also pBis, pOx
[~, ~, pBisModel, pOxModel] = likelihoodGW(p, obsBis, obsOx, EBis, EOx);


end

