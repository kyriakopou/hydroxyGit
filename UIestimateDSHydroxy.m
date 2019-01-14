%parameter estimation procedure
function [xminm, paramNames, paramsOut, sigma, logLik, pData, pModel, pAllStates, ...
    Cov, CovPlot, apprCovFlag, knots, maxDegree, ubR, model] = UIestimateDSHydroxy(dataPoints, obsMatrix, ... 
    errorsMatrix, process, modelFileName)

%turn off annoying warnings
warning('off', 'all');
warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:table:ModifiedVarnames');

global Hessian
%readData from the files Bis.txt and Ox.txt 
%obsMatrix is a matrix of size (numOfDays, numOfObsStates, numOfExperiments)
%errorsMatrix is a matrix of size (numOfDays, size([c, d, e, f], 2), numOfExperiments)
%the size of the matrix is s.t. both dimensions are equal

numOfDays = size(dataPoints, 1);
                
%number of observable states
numObsStates = size(obsMatrix, 2);

%numOfStates depending on what experim. we did
numOfExperiments = sum(any(any(obsMatrix)));
%throw out experiments with all zero observations from obsMatrix
indToRem = find(sum(any(obsMatrix), 2) == 0);
obsMatrix(:,:,indToRem) = [];
errorsMatrix(:,:,indToRem) = [];

%numOfStates depending on what experim. we did
numOfHiddenStates = getNumberOfStatesTemp(numOfExperiments);

%get model's params and construct the matrices of the model and their derivatives if they are not
%constructed yet
[model, paramNames, maxDegree, knots, matrixFileNames] = produceMatricesFiles(modelFileName);
if strcmp(model, 'ctmc')
    matrixFileNames{end} = produceQMatrixFile(maxDegree, knots);
end
numOfParams = length(paramNames);

%function handle to either DSHydroxyEmbryo_ctmc or DSHydroxyEmbryo_dtmc
DSHydroxyEmbryo = str2func(strcat('DSHydroxyEmbryo_', model));


%----TO BE DONE----BASED ON THE PARAMNAMES AND LENGTH OF PARAMS WE HAVE TO CALL HERE CTMC
%SCRIPT TO GET Q AND ITS DERIVATIVES

%initialize output arguments in case we exit earlier (1 timePoint case)
xminm = nan(1, numOfParams);
paramsOut = nan(1, numOfParams);
sigma = nan(1, numOfParams);
fmin = nan;
Cov = nan(numOfParams, numOfParams);
CovPlot = nan(numOfParams, numOfParams);
apprCovFlag = 0;
EMatrix = zeros(numOfHiddenStates, numObsStates, numOfDays, numOfExperiments);
% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
for t=1:numOfDays
    EMatrix(:,:,t,1) = myConvErrorsBisAct(1-errorsMatrix(t,:,1));
    if numOfExperiments > 1 
        EMatrix(:,:,t,2) = myConvErrorsOxAct(1-errorsMatrix(t,:,2));
    end    
    if numOfExperiments > 2    
        EMatrix(:,:,t,3) = myConvErrorsMSIAct(1-errorsMatrix(t,:,3));
    end    
end

%estimate the initial distribution using MLE
%compute initial probability distribution p0
%reshape s.t. the third dimension goes over experiments
% EMatrix = cat(4, EBis, EOx, EMAB);
E = EMatrix(:,:,1,:);
EMatrixForNullTime = reshape(E, [size(E, 1), size(E,2), size(E, 4), size(E, 3)]);
[p0, ~, ~] = maxLikelihood(obsMatrix(1,:,:), EMatrixForNullTime);

lastDay = dataPoints(numOfDays);

%ATTENTION: WHENEVER WE CHANGE THE MAXDEGREE WE HAVE TO RUN AGAIN CTMC AND PRODUCE MATRICES FUNCTIONS 
% maxDegree = 2;
% knots = [7 9];
numOfKnots = size(knots, 2);
numOfCoef = maxDegree + 1 + numOfKnots;
numOfEff = (numOfParams-1) / numOfCoef;

ACond = zeros(2, maxDegree+1);
for i=1:numOfCoef
    if i<= maxDegree+1
        ACond(1,i) = lastDay^(i-1);
        ACond(2,i) = -lastDay^(i-1);
    else    
        ACond(1,i) = (lastDay - knots(i-(maxDegree+1)))^(maxDegree);
        ACond(2,i) = -(lastDay - knots(i-(maxDegree+1)))^(maxDegree);
    end
end

%upper bound for the efficiencies
if strcmp(model, 'dtmc')
    ubR = 1;
else
    %give CTMC upper bound
    ubR = 1;
end

    %0 < = x(1) + t_max*x(2) (+ t_max*x(3)^2) + ... <= 1
A = [ACond, zeros(2, numOfParams - 1*(numOfCoef));
     zeros(2, 1*(numOfCoef)), ACond, zeros(2, numOfParams - 2*(numOfCoef));
     zeros(2, 2*(numOfCoef)), ACond, zeros(2, numOfParams - 3*(numOfCoef));
     zeros(2, 3*(numOfCoef)), ACond, zeros(2, numOfParams - 4*(numOfCoef));
     zeros(2, 4*(numOfCoef)), ACond, zeros(2, numOfParams - 5*(numOfCoef))];

     
b = repmat([ubR; 0], [numOfEff, 1]);

%define the parameter bounds (Initialize both to the same value to fix a parameter)

if (numOfDays >=2 )  

    %step wise function bounds (THIS SHOULD BE TAKEN SYMBOLICALLY)   
    [ubVec, lbVec] = boundsForEachCoef(ubR, lastDay, maxDegree, knots);
    %for the hybrid model ubR for maintenance is 1 and for the rest can be everything
    ubR_maint = 1;   
    [ubVec_maint, lbVec_maint] = boundsForEachCoef(ubR_maint, lastDay, maxDegree, knots);
    
    if strcmp(process{2}, 'rep') 
        switch numOfExperiments
            case 3 
                lb = [repmat(lbVec, [1, numOfEff]), 0];
                ub = [repmat(ubVec, [1, numOfEff]), 1];
            case 2
                lb = [repmat(lbVec, [1, numOfEff-2]), zeros(1, numOfParams - (numOfEff-2)*numOfCoef)];
                ub = [repmat(ubVec, [1, numOfEff-2]), zeros(1, numOfParams - (numOfEff-2)*numOfCoef)];
            otherwise
                lb = [repmat(lbVec, [1, numOfEff-2]), zeros(1, numOfParams - (numOfEff-3)*numOfCoef)];
                ub = [repmat(ubVec, [1, numOfEff-2]), zeros(1, numOfParams - (numOfEff-3)*numOfCoef)];   
        end
    else
        %if no-rep data (WORKS ONLY FOR CONSISTENT NO-REP) only maintenance 0
        lb = [repmat(lbVec, [1, numOfEff]), 0];
        ub = [repmat(ubVec, [1, numOfEff]), 1];
        lbVec_maint = zeros(1:numOfCoef);
        ubVec_maint = zeros(1:numOfCoef);  
    end
    
    lb(1:numOfCoef) = lbVec_maint;
    ub(1:numOfCoef) = ubVec_maint;
    
    %put active demethylation to 0
    lb(numOfCoef*(numOfEff-1)+1:end-1) = zeros(1, numOfCoef);
    ub(numOfCoef*(numOfEff-1)+1:end-1) = zeros(1, numOfCoef);  
    
%if only one time point return the MLE estimators for p0
else
    pAllStates = p0;
    return;
  
end



% Initial guess for parameter vector
x0 = (lb + ub) / 2;      
%incorporate pi0 in the unknown parameters
%constraints for pi0 and modify approprietly the bounds and constraints
% Aeq = [zeros(1, numOfParams), ones(1, numOfHiddenStates)];
% beq = 1;
% lb = [lb, zeros(1, numOfHiddenStates)];
% ub = [ub, ones(1, numOfHiddenStates)];
% x0 = [x0,p0];
% A = [A, zeros(size(A, 1), numOfHiddenStates)];


%for Tet K0 experiments
% lb = [0, -1/lastDay, 0, -1/lastDay, 0, 0, 0, 0, 0.5];
% ub = [1, 1/lastDay, 1, 1/lastDay, 0, 0, 0, 0, 1];


options = optimoptions(@fmincon, 'Display', 'iter-detailed', 'Algorithm', 'interior-point', 'GradObj', 'on', ...
    'MaxFunEvals', 1000, 'AlwaysHonorConstraints', 'bounds', 'StepTolerance', 1e-6,...
    'FunctionTolerance', 1e-6, 'SpecifyConstraintGradient', true, 'FinDiffType', 'central', 'DerivativeCheck', 'off');

% 'FiniteDifferenceStepSize', 1e-8
%reproduce results for the derivative check
% rng(0,'twister'); 

%AlwaysHonorConstraints == 'bounds' then the optimum changes a bit
obj_func = @(x)DSHydroxyEmbryo(x, maxDegree, knots, p0, dataPoints, obsMatrix, EMatrix, ... 
    process, 1, matrixFileNames);
nonlcon_const_func = @(x)ellipseparabola(x, maxDegree, knots, ubR, dataPoints(end));

if maxDegree > 1
    %create the optim problem with the non-linear constraints
    problem = createOptimProblem('fmincon', 'x0', x0, 'objective', obj_func, 'Aineq', A, 'bineq', b, ...
        'lb', lb, 'ub', ub, 'nonlcon', nonlcon_const_func, 'options', options);
else
    %create the optim problem without non-linear constraints
    problem = createOptimProblem('fmincon', 'x0', x0, 'objective', obj_func, 'Aineq', A, 'bineq', b, ...
        'lb', lb, 'ub', ub, 'options', options);
end

%use only fmincon
%[x fval] = fmincon(obj_func, x0, A, b, [], [], lb, ub, [], options);
%use many cores in parallel (it does not give speedup for me)
% currentPool = gcp;
% if(currentPool.Connected == false)
%     parpool(4);
% end

%-----RUN MULTISTART---------
%'StartPointsToRun', 'all' will run exactly k solvers (some of them won't converge)
%'StartPointsToRun', 'bounds-ineqs' will run only solvers that will converge
numOfStartPoints = 1;

globMethod = 'multiStart';

if strcmp(globMethod, 'multiStart')
    %     ms = MultiStart('Display', 'iter', 'StartPointsToRun', 'bounds', 'TolX', 1e-6, 'MaxTime', 1000);
    %     [xminm, fmin, ~, ~, ~] = run(ms, problem, numOfStartPoints);
    % generate custom initial points for multistart that are within bounds
    % but also satisfy the problem linear constraints
    pts = zeros(numOfStartPoints, numOfParams);
    k=1;
    while k<= numOfStartPoints
        ptsTemp = unifrnd(lb, ub);
        if A*ptsTemp' <= b
            pts(k,:) = ptsTemp;
            k = k+1;
        end    
    end
    stpoints = CustomStartPointSet(pts);

    % 'TolX', 1e-10, 'TolFun', 1e-6,
    ms = MultiStart( 'Display', 'iter');
    [xminm, fmin, ~, ~, ~] = run(ms, problem, stpoints);

elseif strcmp(globMethod, 'globSearch')
%-----RUN GLOBALSEARCH------- 
    gs = GlobalSearch('Display', 'iter', 'StartPointsToRun', 'bounds-ineqs', 'MaxTime', 1000, ...
        'NumStageOnePoints', numOfStartPoints, 'NumTrialPoints', 100);
    [xminm, fmin, ~, ~, ~] = run(gs, problem);
end


%call DSHydroxy once more to get the hessian and the prediction of the
%observable and hidden states in case u have computed 2nd derivatives
[~, ~, pData, pModel, pAllStates] = DSHydroxyEmbryo(xminm, maxDegree, knots, p0, dataPoints, ...
    obsMatrix, EMatrix, process, 2, matrixFileNames);
hessian = Hessian;

del = find((ub == lb));
keep = find((ub ~= lb));
hessian(del,:) = [];
hessian(:,del) = [];
paramsOut = xminm;
paramsOut(del) = [];
paramNames(del) = [];

%square root of the observed Fisher Information (i.e. hessian) is an
%approximate lower bound for the standard deviations of the observations
%we dont get -hessian because we compute -logLik
Cov = inv(hessian);

positiveDef = all(eig(hessian) >= 0);
% if (~positiveDef)
%     apprCovFlag = 1;
%     fprintf('Warning! stds are approximations after computing \nthe nearest symmetric-positive definite Cov matrix \n');
% 
%     Cov = nearestSPD(Cov);
% 
% end


%GET THE MATLAB HESSIAN IF SO FAR COV HAS INF OR NAN VALUES
if any(isnan(Cov)) | any(isinf(Cov))
    % %call fmincon again to compute and get the hessian matrix returned by matlab
    x = 0.99 * xminm;
    % %we don't put the constraints here to give more freedom to the parameters
    % %and the model to fit without huge stds (maybe use better fminunc)
    %[xminMatlab, ~, ~, ~, ~, ~, hessianMatlab] = fmincon(obj_func, x, [], [], [], [], lb, ub, [], options);
    %here we impose the constraints
    if maxDegree > 1
        [xminMatlab, ~, ~, ~, ~, ~, hessianMatlab] = fmincon(obj_func, x, A, b, [], [], lb, ub, ...
            nonlcon_const_func, options);
    else
        [xminMatlab, ~, ~, ~, ~, ~, hessianMatlab] = fmincon(obj_func, x, A, b, [], [], lb, ub, [], options);
    end
    
    %throw out of the hessian the nuisance parameters 
    hessianMatlab(del,:) = [];
    hessianMatlab(:,del) = [];

    %square root of the observed Fisher Information (i.e. hessian) is an
    %approximate lower bound for the standard deviations of the observations
    Cov = inv(hessianMatlab);

end

%get sigma from covariance
sigma = (sqrt(diag(Cov)))';    

% create CovPlot putting zero in entries of fixed parameters
% and the Cov entry in entries of non-fixed params
CovPlot = zeros(numOfParams, numOfParams);
for i=1:size(Cov,1)
    for j=1:size(Cov,2)
        CovPlot(keep(i), keep(j)) = Cov(i,j); 
    end    
end    


logLik = -fmin;


%%Likelihood Test for verification
% likStat = zeros(1, 3); 
% pLikTest = zeros(1, 3);
% for testParam = 2:2:6
%     lbNew = lb;
%     ubNew = ub;
%     %define the parameter bounds (Initialize both to the same value to fix a parameter)
%     lbNew(testParam) = 0;
%     ubNew(testParam) = 0;
%     problem = createOptimProblem('fmincon', 'x0', x0, 'objective', obj_func, 'Aineq', A, 'bineq', b, 'lb', lbNew, ...
%   'ub', ubNew,'options', options);
%     ms = MultiStart('StartPointsToRun', 'bounds', 'TolX', 1e-12, 'TolFun', 1e-12);
%     [xminm, fminNull, ~, ~, ~] = run(ms, problem, 10);
%     likStat(testParam / 2) = 2*(fminNull - fmin);
%     pLikTest(testParam / 2) = lratiotest(fmin, fminNull, 1);
% end
% fprintf('The Likel. statistics are: ');
% disp(likStat);







end
