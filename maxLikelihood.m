function [p, pData, pModel] = maxLikelihood(obsMatrix, errorMatrix)
%OXBISOLVE Solve equation system to get initial distribution by estimating
%the maximum likelihood estimator for the initial probability distribution

%numOfStates depending on what experim. we did
numOfExperiments = size(obsMatrix, 3);
numOfStates = getNumberOfStatesTemp(numOfExperiments);

numObsStates = size(obsMatrix(:,:,1), 2);
numOfDays = size(obsMatrix, 1);
%for BS,oxBS,msiBS get the frequencies of the data
%i goes over experiments
pData = zeros(numOfDays, numObsStates, numOfExperiments);
for i=1:numOfExperiments
    pData(:,:,i) = obsMatrix(:,:,i) ./ repmat(nansum(obsMatrix(:,:,i), 2), 1, numObsStates);
end


% Initial guess for probability distribution
% p0 = [0.1, 0.1, zeros(1, numOfStates-2)];
%uniform initial
Aeq = [ones(1,numOfStates); zeros(numOfStates-1, numOfStates)];
beq = [1; zeros(numOfStates-1,1)];

%+mabBS
if numOfExperiments > 2
    p0 = ones(1, numOfStates) * 1/numOfStates; 
    %define the parameter bounds
    lb = zeros(1, numOfStates);
    ub = ones(1, numOfStates);
%+OXbs    
elseif numOfExperiments > 1
    %define the parameter bounds
    lb = zeros(1, numOfStates);
    ub = [1 1 1 0 ones(1,3) 0 ones(1,3) zeros(1, 5)];
    p0 = ub * 1/sum(ub);
    
%only BS
else
    %define the parameter bounds
    lb = zeros(1, 16);
    ub = [1 1 0 0 1 1, zeros(1, 10)];
    p0 = ub * 1/(sum(ub)); 

end


%bounds for david's states


%FOR PARALLEL EXECUTION
% currentPool = gcp;
% if(currentPool.Connected == false)
%     parpool(4);
% end

options = optimset('Algorithm', 'interior-point', 'GradObj', 'on', 'Display', 'off', ...
    'DerivativeCheck', 'off', 'MaxFunEvals', 100);
obj_func = @(p)likelihood(p, obsMatrix, errorMatrix);
problem = createOptimProblem('fmincon', 'x0', p0, 'objective', obj_func, 'Aeq', Aeq, 'beq', beq, ...
    'lb', lb, 'ub', ub, 'options', options);
ms = MultiStart('StartPointsToRun', 'bounds-ineqs', 'Display', 'off', 'MaxTime', 20);
p = run(ms, problem, 5);


%get also prediction of the model
[~, ~, pModel] = likelihood(p, obsMatrix, errorMatrix);


end

