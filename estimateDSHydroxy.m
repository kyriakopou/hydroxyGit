%parameter estimation procedure
%CALL THIS FUNCTION INSTEAD OF GUI WHEN YOU WANT RESULTS FOR MANY REGIONS
%(either singleCpGs or whole regions) -- TO STORE INTO FIXED FOLDER etc.
function [xminm, sigma, fmin] = estimateDSHydroxy(dataPath, dataType)


global hessian
%readData from the files Bis.txt and Ox.txt 
%errorsBis is a matrix of 3 columns (c, d, e) 
%errorsOx is a matrix of 3 columns (c, d, f) 
[dataPointsName, dataPoints, dataLabels, obsBis, obsOx, errorsBis, errorsOx, process] = readData(dataPath);

%keep the region name as a string
slashInd  = strfind(dataPath, '/');
regionName = dataPath(slashInd(end)+1:end);

numOfDays = size(dataPoints, 1);
numObsStates = size(obsBis, 2);
numOfStates = 9;
numOfParams = 9;

% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
EBis = zeros(numOfStates, numObsStates, numOfDays);
EOx = zeros(numOfStates, numObsStates, numOfDays);
for t=1:numOfDays
    EBis(:,:,t) = convErrorBis(errorsBis(t,:));
    EOx(:,:,t) = convErrorOx(errorsOx(t,:));   
end

lastDay = dataPoints(numOfDays);


%cell array containing paramNames
paramNames = {'b0_m', 'b1_m', 'b0_d', 'b1_d', 'b0_e', 'b1_e', 'b0_dm', 'b1_dm', 'p'};
% Initial guess for parameter vector
x0 = [0.5, -0.5/lastDay, 0.1, 0.1/lastDay, 0.1, 0.1/lastDay, 0, 0, 0.7];

%estimate the initial distribution using lsqnonneg
% pBis0 = obsBis(1,:) / sum(obsBis(1,:), 2);
% pOx0 = obsOx(1,:) / sum(obsOx(1,:), 2);
% p0 = oxbisolve(pBis0, pOx0, EBis(:,:,1), EOx(:,:,1)); 
% p0 = p0';

%estimate the initial distribution using MLE
%compute initial probability distribution p0
p0 = maxLikelihood(obsBis(1,:), obsOx(1,:), EBis(:,:,1), EOx(:,:,1), dataType);

%constraints for x(1), x(2)
%x(3), x(4), x(5), x(6)
maxDegree = 1;
ACond = zeros(2, maxDegree+1);
for i=1:maxDegree+1 
    ACond(1, i) = lastDay^(i-1);
    ACond(2, i) = -lastDay^(i-1);
end


    %0 < = x(1) + t*x(2) <= 1
A = [ACond, zeros(2, numOfParams-(maxDegree+1));
    %0 < = x(3) + t*x(4) <= 1
     zeros(2, 1*(maxDegree+1)), ACond, zeros(2, numOfParams - 2*(maxDegree+1));
     %0 < = x(5) + t*x(6) <= 1
     zeros(2, 2*(maxDegree+1)), ACond, zeros(2, numOfParams - 3*(maxDegree+1));
     %0 < = x(7) + t*x(8) <= 1
     zeros(2, 3*(maxDegree+1)), ACond, zeros(2, numOfParams - 4*(maxDegree+1))];
     
     
b = repmat([1; 0], 4, 1);

%define the parameter bounds (Initialize both to the same value to fix a parameter)
if (strcmp(dataType, 'oxBS'))
    if (numOfDays > 2)
        %linear maintenance, de novo, hydroxy
        lb = [0, -1/lastDay, 0, -1/lastDay, 0, -1/lastDay, 0, 0, 0];
        ub = [1, 1/lastDay, 1, 1/lastDay, 1, 1/lastDay, 0, 0, 1];
    else
        %constant maintenance, de novo, hydroxy
        lb = [0, 0, 0, 0, 0, 0, 0, 0, 0.5];
        ub = [1, 0, 1, 0, 1, 0, 0, 0, 1];
    end    
else
    if (numOfDays > 2)
        %constant maintenance , linear de novo
        lb = [0, 0, 0, -1/lastDay, 0, 0, 0, 0, 0];
        ub = [1, 0, 1, 1/lastDay, 0, 0, 0, 0, 0];
    else
        %constant maintenance, de novo
        lb = [0, 0, 0, 0, 0, 0, 0, 0, 0];
        ub = [1, 0, 1, 0, 0, 0, 0, 0, 0];
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
    'HessFcn', @hessFcn, 'MaxFunEvals', 20000, 'AlwaysHonorConstraints', 'bounds',  'FinDiffType','central', 'DerivativeCheck', 'off');
obj_func = @(x)DSHydroxy(x, p0, dataPoints, obsBis, obsOx, EBis, EOx, process);
problem = createOptimProblem('fmincon', 'x0', x0, 'objective', obj_func, 'Aineq', A, 'bineq', b, 'lb', lb, 'ub', ub,'options', options);
% ms = MultiStart('StartPointsToRun', 'bounds', 'TolX', 1e-12, 'TolFun', 1e-12);
ms = MultiStart('StartPointsToRun', 'bounds', 'TolX', 1e-10, 'TolFun', 1e-6);
[xminm, fmin, ~, output, manyminsm] = run(ms, problem, 10);


%gs = GlobalSearch('StartPointsToRun', 'bounds-ineqs', 'TolX', 1e-12);
% [xminm,fminm,flagm,outptm,manyminsm] = run(gs, problem);

%call DSHydroxy with the optimal parameter values once more to get the hessian 
%and the model prediction for each data point
[~, ~, pBis, pOx, pBisModel, pOxModel, pAllStates] = DSHydroxy(xminm, p0, dataPoints, obsBis, obsOx, EBis, EOx, process);

% while (numel(xminm) == 0)
%     try
%         [~, ~, pBis, pOx, pBisModel, pOxModel, pAllStates] = DSHydroxyEmbryo(xminm, p0, dataPoints, obsBis, obsOx, EBis, EOx, process);
%     catch ME
%         disp('Restarting Multistart ...');
%         [xminm, fmin, ~, output, manyminsm] = run(ms, problem, 2);
%     end
% end

H = hessian;
%throw out of the hessian the nuisance parameters 
del = find((ub == lb));
keep = find((ub ~= lb));
H(del,:) = [];
H(:,del) = [];
paramsOut = xminm;
paramsOut(del) = [];
paramNames(del) = [];


%call sensitivity analysis function
%sensitivityAnalysis(xminm, p0, dataPoints, obsBis, obsOx, EBis, EOx, process, pAllStates, paramNames, regionName);


%square root of the Fisher Information(of the observations) is an
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

%if previous covariance is not positive semi-definite
%-> has negative diagonal elements
positiveDef = all(eig(H) >= 0);
if (~positiveDef)
    
    fprintf('Warning: stds are computed approximately by fmincon hessian \n because the numerical hessian was not positive semi-definite \n');
    %call fmincon again to compute and return the hessian matrix by matlab
    x = 0.98 * xminm;
    % %we don't put the conditions here to give more freedom to the parameters
    % %and the model to fit without huge stds (maybe use better fminunc)
    [xminMatlab, ~, ~, ~, ~, ~, hessianMatlab] = fmincon(obj_func, x, [], [], [], [], lb, ub, [], options);
    % %throw out of the hessian the nuisance parameters 
    i = find((ub==0 & lb==0) | (ub==1 & lb ==1));
    hessianMatlab(i,:) = [];
    hessianMatlab(:,i) = [];
    CovMatlab = inv(hessianMatlab);
    sigmaMatlab = (sqrt(diag(CovMatlab)))';
    sigma = sigmaMatlab;

end



%display(manyminsm);
% 
% n = size(manyminsm, 2);
% figure('name', 'Estimates');
% for i=1:n;
%   subplot(3,1,1), plot(i, manyminsm(1,i).X(1),'o');
%   hold on;
%   subplot(3,1,1), plot(i, manyminsm(1,i).X(2),'+');
%   hold on;
%   subplot(3,1,2), plot(i, manyminsm(1,i).X(3),'*');
%   hold on;
% end

%   fcompl = @(y) f(y)-f(x);
%   display(f(x));
%   [rnd, eval] = slicesample(x0, 10000, 'pdf', f);%'width', [0.5 0.5 0.5 0.5 0.5 0.5 0.5]);
%   display(rnd);
%   %display(rnd(1,:));
%   C = num2cell(transpose(rnd),1);
%   %display(C);
%   res = cellfun(fcompl,C);
%   display(res);
%   display(f(x));
% %   [binheight,bincenter] = hist(rnd(:,7),50);
% %   h = bar(bincenter,binheight,'hist');
% %   set(h,'facecolor',[0.8 .8 1]);
%   %plotmatrix(exp(rnd(:,1)),exp(rnd(:,7)));
%   scatter(rnd(:,7),res);


%count all local optimum different than the global one
localButNotGlobal = 0;
for i=1:size(manyminsm, 2)
    if(~isequal(round(manyminsm(1,1).X, 2), round(manyminsm(1,i).X, 2)))
        localButNotGlobal = localButNotGlobal + 1;
    end
end


%call the main output function. xminm contains all the parameters and
%CovPlot has zeros in entries of nuisance params
%fig = plotOutput(xminm, p0, dataPoints, obsBis, obsOx, EBis, EOx, CovPlot, dataPointsName, labels, process);

if (size(regionName, 2) < 4) || (~strcmp(regionName(end-3:end-1), 'CpG'))
    fig = plotOutput(xminm, pBis, pOx, pBisModel, pOxModel, pAllStates, CovPlot, dataPointsName, dataPoints, dataLabels, dataType, regionName);
else
    fig = plotOutputSingleCpG(xminm, pBis, pOx, pBisModel, pOxModel, pAllStates, CovPlot, dataPointsName, dataPoints, dataLabels, dataType, regionName);

end    
%fig = plotRates(xminm, dataPoints, CovPlot, regionName);

set(fig, 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape', 'PaperType', 'a4');
picPath = strcat('/Users/kyriakopou/Desktop/', regionName);
print(picPath, '-dpdf', '-r0')
%plotRates(xminm, days, CovPlot, regionName);

%THIS PATH HAS TO BE GIVEN AS AN ARGUMENT BY THE USER!
outFilePath = strcat('/Users/kyriakopou/Desktop/', regionName, '.txt');

fileID = fopen(outFilePath, 'w');

fprintf('\n');
disp('Optimization is done! Check the file  ')
disp(outFilePath);
fprintf('\n');

fprintf(fileID, '%s \n', output.message);
fprintf(fileID, 'The local but not global maxima are: %d \n\n', localButNotGlobal);

%Print the estimated parameters the stds
%and relative half widths
fprintf(fileID, 'The parameters of the model are: \n\n');
formatSpec = '%s = %4.4f ? (%4.4f)\n';
for i = 1:size(paramsOut, 2)
    %fprintf(fileID, formatSpec, paramNames{1,:});
    %formatSpec = '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f\n';
    fprintf(fileID, formatSpec, char(paramNames(i)), paramsOut(i), sigma(i));
end    
fprintf(fileID, '\n');



% fprintf('The relative half widths are: \n');
% relHalfWidth = abs(1.96*sigma ./ params);
% display(relHalfWidth);
% fprintf('\n');

lambdaCoef(paramsOut, Cov, fileID)
fprintf(fileID, '-------------------------------------\n');

if (numOfDays > 2)
    if (strcmp(dataType, 'oxBS'))

        %define the g function and the Jacobian for the rate params of interest (not-nuisance)
        Wald = zeros(1, floor(size(paramsOut, 2) / 2));
        for testParam = 2:2:size(paramsOut, 2)
            g = paramsOut(testParam);
            Jacob = zeros(1, size(paramsOut,2));
            Jacob(testParam) = 1;
            Wald(testParam / 2) = g * inv(Jacob * Cov * Jacob')*g';
            %convert cell array position to string (char array)
            paramName = char(paramNames(testParam));
            fprintf(fileID, 'The Wald statistic for the %s parameter is: %d \n', paramName, Wald(testParam / 2));
        end

        %wald test for non-recognition prob. if it is in the keeped params
        if (strcmp(paramNames{end}, 'p') )
            gP = paramsOut(end) - 1;
            JacobP = zeros(1, size(paramsOut, 2));
            JacobP(7) = 1;
            WaldP = gP * inv(JacobP*Cov*JacobP')*gP';
            fprintf(fileID, 'The Wald statistic for the p parameter is: %d\n', WaldP);
        end

        %wald test for total methylation lambda
        gLambda = [paramsOut(2) + paramsOut(4) - paramsOut(2)*paramsOut(3) - paramsOut(1)*paramsOut(4);
                    paramsOut(2)*paramsOut(4)];
        JacobLambda = [-paramsOut(4), 1-paramsOut(3), -paramsOut(2), 1-paramsOut(1), zeros(1, size(paramsOut, 2)-4);
                       0, paramsOut(4), 0, paramsOut(2), zeros(1, size(paramsOut, 2) - 4)];
        WaldLambda = gLambda' * inv(JacobLambda*Cov*JacobLambda')*gLambda;
        fprintf(fileID, 'The Wald statistic for total methylation is: %d\n', WaldLambda);

    else

        %define the g function and the Jacobian for the rate params of
        %interest (here it is only b1_denovo)
        g = paramsOut(3);
        Jacob = zeros(1, size(paramsOut,2));
        Jacob(3) = 1;
        Wald = g * inv(Jacob * Cov * Jacob')*g';
        %convert cell array position to string (char array)
        paramName = char(paramNames(3));
        fprintf(fileID, 'The Wald statistic for the %s parameter is: %d \n', paramName, Wald);

        %wald test for total methylation lambda
        gLambda = paramsOut(3) - paramsOut(1)*paramsOut(3);

        JacobLambda = [-paramsOut(3), 0, 1-paramsOut(1)];
        WaldLambda = gLambda' * inv(JacobLambda*Cov*JacobLambda')*gLambda;
        fprintf(fileID, 'The Wald statistic for total methylation is: %d\n', WaldLambda);

    end
end
    
    fprintf(fileID, '---------------------------------------------------------------\n');
    
    C = num2cell(pBis);
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';
    
    fprintf(fileID, 'The data distribution of BS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec,C{row,:});
    end
    fprintf(fileID, '\n');
    
    C = num2cell(pBisModel);
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';
    
    fprintf(fileID, 'The predicted distribution of BS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec,C{row,:});
    end
    
    fprintf(fileID, '---------------------------------------------------------------\n');
    
    C = num2cell(pOx);
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';
    
    fprintf(fileID, 'The data distribution of oxBS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec,C{row,:});
    end
    fprintf(fileID, '\n');
    
    C = num2cell(pOxModel);
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';
    
    fprintf(fileID, 'The predicted distribution of oxBS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec,C{row,:});
    end
    
    fprintf(fileID, '---------------------------------------------------------------\n');
    
    C = num2cell(pAllStates);
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n';
    
    fprintf(fileID, 'The predicted distribution of the hidden states \n\n');
    fprintf(fileID, 'data uu um mu uh hu hm mh mm hh \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec,C{row,:});
    end
    fprintf(fileID, '\n');
    
    


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
