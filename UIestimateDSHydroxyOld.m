%parameter estimation procedure
function [xminm, paramNames, paramsOut, sigma, fmin, pBis, pOx, pBisModel, pOxModel, pAllStates, ...
    Cov, CovPlot, apprCovFlag] = UIestimateDSHydroxy(dataPoints, obsMatrix, errorsMatrix, process, regionName, dataPointsName, dataLabels)

%turn off annoying warnings
warning('off', 'all');
warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:table:ModifiedVarnames');

global hessian
%readData from the files Bis.txt and Ox.txt 
%errorsBis is a matrix of 3 columns (c, d, e) 
%errorsOx is a matrix of 3 columns (c, d, f) 

% BisFile = strcat(regionName, '_BS.txt');
% OxFile = strcat(regionName, '_oxBS.txt');
% obsBis = samplingData(BisFile);
% obsOx = samplingData(OxFile);

numOfDays = size(dataPoints, 1);

%AFP
% obsMatrix(:,:,3) = [0, 43, 22, 1122; 35, 118, 173, 1862; 0, 81, 73, 1274; 8, 92, 98, 1374]; 
% errorsVectors(:,:,3) = [0.0101967880085653, 0.0728051391862955, 0.0235546038543897;
%                         0.0101967880085653, 0.0728051391862955, 0.0235546038543897
%                         0.00806759868421048, 0.125, 0.121710526315789;
%                         0.00817643504531718, 0.0211480362537764, 0.11178247734139];
%IAP
obsMatrix(:,:,3) = [32, 277, 279, 3262;18, 159, 179, 2137; 5, 89, 89, 1036;17, 145, 154, 1956];
errorsVectors(:,:,3) =  [0.0114167528438469, 0.047569803516029, 0.0496380558428128;
                         0.0118759936406996, 0.0763116057233704, 0.0508744038155803;
                         0.0135830618892508, 0.0651465798045603, 0.0586319218241042;
                         0.0125989672977624, 0.044750430292599, 0.0464716006884682];
% %dmr1
% obsMatrix(:,:,3) = [11, 63, 215, 2470; 5, 98, 54, 2180; 1, 98, 37, 1833;2, 127, 92, 1990];
% errorsVectors(:,:,3) = [0.00905900900900902, 0.0630630630630631, 0; 
%                         0.0084913978494624, 0.0161290322580645, 0;
%                         0.00841455696202531, 0.145569620253165, 0;
%                         0.00688114285714281, 0.0514285714285714, 0];
% %L1
% obsMatrix(:,:,3) = [396, 4290, 4544, 34841; 1015, 6173, 5734, 46743; 857, 7133, 7159, 55738; 372, 2447, 2060, 18781];
% errorsVectors(:,:,3) = [0.0100844898400346, 0.122243839169909, 0.10516645049719;
%                         0.0104716027321102, 0.112064172821857, 0.0972917163052974;
%                         0.00993475923328657, 0.129499766245909, 0.102517865491218;
%                         0.0113026606425702, 0.0873493975903614, 0.112248995983936];
% %mSat
% obsMatrix(:,:,3) = [35, 235, 210, 1967; 30, 174, 153, 1551; 14, 135, 136, 1235; 24, 189, 162, 1752];
% errorsVectors(:,:,3) = [0.00947033898305083, 0.00529661016949153, 0.230932203389831;
%                         0.0114864864864865, 0.00405405405405405, 0.218918918918919;
%                         0.00993197278911562, 0.00510204081632653, 0.215986394557823;
%                         0.00860943168077388, 0.00846432889963724, 0.216444981862152];
% %muERVL
% obsMatrix(:,:,3) = [155, 1342, 1312, 13091; 147, 1577, 1919, 14827; 122, 1398, 1473, 12485;179, 1303, 1284, 12636];
% errorsVectors(:,:,3) = [0.00946718631897203, 0.00245653817082389, 0.24244142101285;
%                         0.00954133903602561, 0.0021385096232933, 0.22174699786149;
%                         0.00888065229940349, 0.00384837406195882, 0.280931306522994;
%                         0.00827074644318848, 0.00331319430910154, 0.244591697524849];
%DMR2
% obsMatrix(:,:,3) = [5, 409, 210, 104870; 
%                     3, 405, 270, 87281; 
%                     137, 8172, 4486, 58559; 
%                     35, 4286, 5017, 118290; 
%                     30, 1233, 5715, 45550; 
%                     123, 4106, 3992, 63957; 
%                     20, 2289, 1563, 86749];


%number of observable states
numObsStates = size(obsMatrix, 2);


%numOfStates depending on what experim. we did
numOfExperiments = size(obsMatrix, 3);
numOfStates = getNumberOfStatesTemp(numOfExperiments);


errorsBis = [errorsMatrix(:,:,1), zeros(numOfDays, 1)];
errorsOx = [errorsMatrix(:,:,2), zeros(numOfDays, 1)];
errorsMAB = [errorsMatrix(:,:,3), repmat(0.065, [numOfDays, 1])];

%cell array containing paramNames
paramNames = {'b0_m', 'b1_m', 'b0_d', 'b1_d', 'b0_e', 'b1_e', 'b0_dm', 'b1_dm', 'p'};
numOfParams = length(paramNames);

%initialize output arguments in case we exit earlier (1 timePoint case)
xminm = nan(1, numOfParams);
paramsOut = nan(1, numOfParams);
sigma = nan(1, numOfParams);
fmin = nan;
Cov = nan(7, 7);
CovPlot = nan(numOfParams, numOfParams);
apprCovFlag = 0;

% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
EBis = zeros(numOfStates, numObsStates, numOfDays);
EOx = zeros(numOfStates, numObsStates, numOfDays);
EMAB = zeros(numOfStates, numObsStates, numOfDays);
for t=1:numOfDays
    EBis(:,:,t) = myConvErrorsBisAct(1-errorsBis(t,:));
    EOx(:,:,t) = myConvErrorsOxAct(1-errorsOx(t,:));
    EMAB(:,:,t) = myConvErrorsMSIAct(1-errorsMAB(t,:));
end


%estimate the initial distribution using MLE
%compute initial probability distribution p0
%reshape s.t. the third dimension goes over experiments
EMatrix = cat(4, EBis, EOx, EMAB);
E = EMatrix(:,:,1,:);
EMatrixForNullTime = reshape(E, [size(E, 1), size(E,2), size(E, 4), size(E, 3)]);
[p0, ~, ~] = maxLikelihood(obsMatrix(1,:,:), EMatrixForNullTime);

lastDay = dataPoints(numOfDays);



%constraints for x(1), x(2)
%x(3), x(4), x(5), x(6)
maxDegree = 1;
ACond = zeros(2, maxDegree+1);
for i=1:maxDegree+1 
    ACond(1, i) = lastDay^(i-1);
    ACond(2, i) = -lastDay^(i-1);
end


    %0 < = x(1) + t_max*x(2) <= 1
A = [ACond, zeros(2, numOfParams-(maxDegree+1));
    %0 < = x(3) + t_max*x(4) <= 1
     zeros(2, 1*(maxDegree+1)), ACond, zeros(2, numOfParams - 2*(maxDegree+1));
     %0 < = x(5) + t_max*x(6) <= 1
     zeros(2, 2*(maxDegree+1)), ACond, zeros(2, numOfParams - 3*(maxDegree+1));
     %0 < = x(7) + t_max*x(8) <= 1
     zeros(2, 3*(maxDegree+1)), ACond, zeros(2, numOfParams - 4*(maxDegree+1))];
     
     
b = repmat([1; 0], 4, 1);

%define the parameter bounds (Initialize both to the same value to fix a parameter)
if numOfExperiments == 3
    if (numOfDays > 2)
        %linear maintenance, de novo, hydroxy
%         lb = [0, -1/lastDay, 0, -1/lastDay, 0, -1/lastDay, 0, -1/lastDay, 0.5];
%         ub = [1,  1/lastDay, 1,  1/lastDay, 1,  1/lastDay, 1,  1/lastDay, 1];
%         % Initial guess for parameter vector
%         x0 = [0.2, -0.2/lastDay, 0.1, 0.1/lastDay, 0.1, 0.1/lastDay, 0.1, 0.1/lastDay, 0.8];
        %linear maintenance, de novo, hydroxy
        lb = [repmat([0, -1/lastDay], [1,(numOfParams-1)/2]), 0];
        ub = [repmat([1,  1/lastDay], [1,(numOfParams-1)/2]), 1];
        % Initial guess for parameter vector
        x0 = [repmat([0.2, 0.2/lastDay], [1,(numOfParams-1)/2]), 0.8];
        
    elseif numOfDays == 2
        %constant maintenance, de novo, hydroxy, active demethylation
        lb = zeros(1, numOfParams);
        ub = [repmat([1,0], [1,(numOfParams-1)/2]), 1];
        % Initial guess for parameter vector
        x0 = [repmat([0.2,0], [1,(numOfParams-1)/2]), 0.8];
    else
        %if only one time point return the MLE estimators for  p0
        pAllStates = p0;
        return;
    end
elseif numOfExperiments == 2    
    if (numOfDays > 2)
        %linear maintenance, de novo, hydroxy, 0 active demethylation
        lb = [repmat([0, -1/lastDay], [1, 3]), zeros(1, numOfParams-6)];
        ub = [repmat([1,  1/lastDay], [1, 3]), zeros(1, numOfParams-7), 1];
        % Initial guess for parameter vector
        x0 = [repmat([0.2, 0.2/lastDay], [1, 3]), zeros(1, numOfParams-7), 0.8];
    elseif numOfDays == 2
        %constant maintenance, de novo, hydroxy, 0 active demethylation 
        lb = zeros(1, numOfParams);
        ub = [repmat([1,0], [1,3]), zeros(1, numOfParams-7), 1];
        % Initial guess for parameter vector
        x0 = [repmat([0.2,0], [1,3]), zeros(1, numOfParams-7), 0.8];
    else
        %if only one time point return the MLE estimators for p0
        pAllStates = p0;
        return;
    end    
else
    if (numOfDays > 2)
        %constant maintenance , linear de novo
%         lb = [0, 0, 0, -1/lastDay, 0, 0, 0, 0, 0];
%         ub = [1, 0, 1, 1/lastDay, 0, 0, 0, 0, 0];       
        %linear maintenance , linear de novo, others 0
        lb = [repmat([0,-1/lastDay], [1,2]), zeros(1, numOfParams-4)];
        ub = [repmat([1,1/lastDay], [1,2]), zeros(1, numOfParams-4)];
        % Initial guess for parameter vector
        x0 = [repmat([0.2,0.2/lastDay], [1, 2]), zeros(1, numOfParams-4)];
    elseif numOfDays == 2
        %constant maintenance, de novo, others 0
        lb = zeros(1, numOfParams);
        ub = [repmat([1,0], [1,2]), zeros(1, numOfParams-4)];
        % Initial guess for parameter vector
        x0 = [repmat([0.1,0], [1,2]), zeros(1, numOfParams-4)];
    else         
        %if only one time point return the MLE estimators for  p0
        pAllStates = p0;
        return;
    end   
end

    
%incorporate pi0 in the unknown parameters
%constraints for pi0
% Aeq = [zeros(1,9), ones(1,9)];
% beq = 1;
% lb = [0, -1/lastDay, 0, -1/lastDay, 0, -1/lastDay, 0, 0, 0.5, zeros(1,9)];
% ub = [1, 1/lastDay, 1, 1/lastDay, 1, 1/lastDay, 0, 0, 1, ones(1,9)];
% x0 = [x0,p0];

%for Tet K0 experiments
% lb = [0, -1/lastDay, 0, -1/lastDay, 0, 0, 0, 0, 0.5];
% ub = [1, 1/lastDay, 1, 1/lastDay, 0, 0, 0, 0, 1];


options = optimset('Display', 'off', 'Algorithm', 'interior-point', 'GradObj', 'on', 'Hessian','user-supplied', ...
    'HessFcn', @hessFcn, 'MaxFunEvals', 100, 'AlwaysHonorConstraints', 'bounds',  'FinDiffType','central', 'DerivativeCheck', 'off');
%if AlwaysHonorConstraints == 'bounds' then the optimum changes a bit
obj_func = @(x)DSHydroxyEmbryo(x, p0, dataPoints, obsMatrix, EMatrix, process, 1);
problem = createOptimProblem('fmincon', 'x0', x0, 'objective', obj_func, 'Aineq', A, 'bineq', b, 'lb', lb, 'ub', ub, 'options', options);

%use only fmincon
%[x fval] = fmincon(obj_func, x0, A, b, [], [], lb, ub, [], options);
%use many cores in parallel (it does not give speedup for me)
% currentPool = gcp;
% if(currentPool.Connected == false)
%     parpool(4);
% end
%'StartPointsToRun', 'all' will run exactly k solvers (some of them won't converge)
%'StartPointsToRun', 'bounds-inew' will run only solvers that will converge
ms = MultiStart('Display', 'iter', 'StartPointsToRun', 'bounds-ineqs', 'TolX', 1e-6, 'MaxTime', 1000);
[xminm, fmin, ~, ~, ~] = run(ms, problem, 20);

%generate custom initial points for multistart that are within bounds
%but also satisfy the problem constraints
% pts = zeros(10, 9);
% k=1;
% while k<=10
%     ptsTemp = unifrnd(lb, ub);
%     if A*ptsTemp' <=b;
%         pts(k,:) = ptsTemp;
%         k = k+1;
%     end    
% end
% stpoints = CustomStartPointSet(pts);
%ms = MultiStart('TolX', 1e-10, 'TolFun', 1e-6, 'Display', 'iter');
% [xminm, fmin, ~, ~, ~] = run(ms, problem, stpoints);

% gs = GlobalSearch('StartPointsToRun', 'bounds-ineqs', 'TolX', 1e-10, 'TolFun', 1e-6, 'Display', 'iter');
% [xminm,fmin, flagm,outptm,manyminsm] = run(gs, problem);

%call DSHydroxy once more to get the hessian and the prediction of the
%observable and hidden states 
[~, ~, pData, pModel, pAllStates] = DSHydroxyEmbryo(xminm, p0, ... 
    dataPoints, obsMatrix, EMatrix, process, 1);

pBis = pData(:,:,1);
pOx = pData(:,:,2);
pBisModel = pData(:,:,1);
pOxModel = pData(:,:,2);

H = hessian;
%throw out of the hessian the nuisance parameters 
del = find((ub == lb));
keep = find((ub ~= lb));
H(del,:) = [];
H(:,del) = [];
paramsOut = xminm;
paramsOut(del) = [];
paramNames(del) = [];

%square root of the observed Fisher Information (i.e. hessian) is an
%approximate lower bound for the standard deviations of the observations
Cov = inv(H);
sigma = (sqrt(diag(Cov)))';

%create CovPlot putting zero in entries of fixed parameters
%and the Cov entry in entries of non-fixed params
CovPlot = zeros(numOfParams, numOfParams);
for i=1:size(Cov,1)
    for j=1:size(Cov,2)
        CovPlot(keep(i), keep(j)) = Cov(i,j); 
    end    
end    

if all(~isnan(Cov)) & all(~isinf(Cov))
    
    %if the hessian is not positive definite -> not invertible
    %(inv function works but gives negative diagonal elements in cov)
    positiveDef = all(eig(H) >= 0);
    if (~positiveDef)
        apprCovFlag = 1;
        fprintf('Warning! stds are approximations after computing \nthe nearest symmetric-positive definite Cov matrix \n');

        pseudoCov = nearestSPD(Cov);
        sigma = sqrt(diag(pseudoCov))';

        %call fmincon again to compute and return the hessian matrix by matlab
    %     x = 0.98 * xminm;
    %     % %we don't put the conditions here to give more freedom to the parameters
    %     % %and the model to fit without huge stds (maybe use better fminunc)
    %     %[xminMatlab, ~, ~, ~, ~, ~, hessianMatlab] = fmincon(obj_func, x, [], [], [], [], lb, ub, [], options);
    %     [xminMatlab, ~, ~, ~, ~, ~, hessianMatlab] = fmincon(obj_func, x0, A, b, [], [], lb, ub, [], options);


        % %throw out of the hessian the nuisance parameters 
    %     i = find((ub==0 & lb==0) | (ub==1 & lb ==1));
    %     hessianMatlab(i,:) = [];
    %     hessianMatlab(:,i) = [];
    %     pseudoInvHessian = pinv(hessianMatlab);
    %     pseudoHessian = nearestSPD(hessianMatlab);

    %     hessianMatUnc(i,:) = [];
    %     hessianMatUnc(:,i) = [];

    %     CovMatlab = inv(hessianMatlab);
    %     CovMatUnc = inv(hessianMatUnc);
    %     pseudoCov = chol(pseudoInvHessian)' * chol(pseudoInvHessian);

    %     sigmaMatlab = (sqrt(diag(CovMatlab)))';    
    %     sigmaMatUnc = (sqrt(diag(CovMatUnc)))';
    %     sigma = sigmaMatlab;


    end

end


% %call the main output function. xminm contains all the parameters and
% %CovPlot has zeros in entries of nuisance params
fig = plotOutput(xminm, pData, pModel, pAllStates, CovPlot, dataPointsName, dataPoints, dataLabels);
%export the figure as pdf in a3 format
set(fig, 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape', 'PaperType', 'a3');
picPath = strcat('/Users/kyriakopou/Desktop/', regionName);
print(fig, picPath, '-dpdf', '-r0')



%verification that wald statistic computation works
%the next does not work for hessian taken by hand because Cov is not symmetric
%probably due to negligible numerical differences
%[h, p, Wald1, crit] = waldtest(g, Jacob, Cov);
%fprintf('The Wald statistic calling matlab function is: %d\n', Wald1);


%%Likelihood Test for verification
% likStat = zeros(1, 3); 
% pLikTest = zeros(1, 3);
% for testParam = 2:2:6
%     lbNew = lb;
%     ubNew = ub;
%     %define the parameter bounds (Initialize both to the same value to fix a parameter)
%     lbNew(testParam) = 0;
%     ubNew(testParam) = 0;
%     problem = createOptimProblem('fmincon', 'x0', x0, 'objective', obj_func, 'Aineq', A, 'bineq', b, 'lb', lbNew, 'ub', ubNew,'options', options);
%     ms = MultiStart('StartPointsToRun', 'bounds', 'TolX', 1e-12, 'TolFun', 1e-12);
%     [xminm, fminNull, ~, ~, ~] = run(ms, problem, 10);
%     likStat(testParam / 2) = 2*(fminNull - fmin);
%     pLikTest(testParam / 2) = lratiotest(fmin, fminNull, 1);
% end
% fprintf('The Likel. statistics are: ');
% disp(likStat);

% fprintf('The Wald statistics are: ');
% disp(Wald);



end
