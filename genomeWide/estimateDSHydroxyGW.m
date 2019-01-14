%parameter estimation procedure
%CALL THIS FUNCTION INSTEAD OF GUI WHEN YOU WANT RESULTS FOR MANY REGIONS
%(either singleCpGs or whole regions) -- TO STORE INTO FIXED FOLDER etc.
function [paramsOut, paramNames, CovOut, pAllStates, KLBis, KLOx, volume] = estimateDSHydroxyGW(chr, binAndCpGs, ... 
    maxObsDaysVector, maxDataLabels, obsBis, obsOx, errorsBis, errorsOx, process, chrDir)

global hessian


%days with non zeros observations
dayInd = find(~all(obsBis == 0, 2));
obsDays = maxObsDaysVector(dayInd);
dataLabels = maxDataLabels(dayInd);

numOfDays = size(obsDays, 1);
numObsStates = size(obsBis, 2);
numOfStates = 9;
numOfParams = 7;
numOfSamplesBS = sum(obsBis, 2);
numOfSamplesoxBS = sum(obsOx, 2);
numOfSamples = sum(numOfSamplesBS + numOfSamplesoxBS);
CovOut = nan(7, 7);
paramsOut = nan(1, 7);
paramNames = nan(1, 7);
pAllStates = nan(3, 9);
pBis = nan(3, 4);
pOx = nan(3, 4);
pBisModel = nan(3, 4);
pOxModel = nan(3, 4);
volume = nan;
SPDFlag = 0;
del = [];
keep = [1:7];

bin = cell2mat(binAndCpGs(1));
cpgArray = cell2mat(binAndCpGs(2:end));
stringCpGs = sprintf('%d_', cpgArray);
regionName = strcat('bin_', int2str(bin));
dataType = 'oxBS';


% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
EBis = zeros(numOfStates, numObsStates, numOfDays);
EOx = zeros(numOfStates, numObsStates, numOfDays);
for t=1:3
    EBis(:,:,t) = convErrorBis(errorsBis(t,:));
    EOx(:,:,t) = convErrorOx(errorsOx(t,:));   
end


%define the parameter bounds (Initialize both to the same value to fix a parameter)
if (numOfDays == 3)
    
    %cell array containing paramNames
    paramNames = {'b0_m', 'b1_m', 'b0_d', 'b1_d', 'b0_e', 'b1_e', 'p'};
    % lastDay = dataPoints(numOfDays);
    lastDay = obsDays(numOfDays);
    % Initial guess for parameter vector
    x0 = [0.5, 0, 0.1, 0.1/lastDay, 0.1, 0.1/lastDay, 0.7];
    %tet_k0
%     x0 = [0.5, -0.5/lastDay, 0.1, 0.1/lastDay, 0, 0, 0, 0, 0];
    %estimate the initial distribution using MLE
    %compute initial probability distribution p0
    p0 = maxLikelihoodGW(obsBis(1,:), obsOx(1,:), EBis(:,:,1), EOx(:,:,1), dataType);

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
         zeros(2, 2*(maxDegree+1)), ACond, zeros(2, numOfParams - 3*(maxDegree+1))];
   

    b = repmat([1; 0], 3, 1);
    
    %linear maintenance, de novo, hydroxy
    lb = [0, -1/lastDay, 0, -1/lastDay, 0, -1/lastDay, 0];
    ub = [1, 1/lastDay, 1, 1/lastDay, 1, 1/lastDay, 1];
    
    %TET-KO, linear maintenance, de novo, no hydroxy
%     lb = [0, -1/lastDay, 0, -1/lastDay, 0, 0, 0, 0, 1];
%     ub = [1, 1/lastDay, 1, 1/lastDay, 0, 0, 0, 0, 1];
    
    options = optimoptions(@fmincon, 'Display', 'off', 'Algorithm', 'interior-point', 'GradObj', 'on', ...
        'Hessian', 'user-supplied', 'HessFcn', @hessFcn, 'MaxFunEvals', 100, 'MaxIterations', 500, ... 
        'OptimalityTolerance', 1e-10, 'FinDiffType', 'central', 'DerivativeCheck', 'off');
% 'AlwaysHonorConstraints', 'bounds',
    obj_func = @(x)DSHydroxyEmbryoGW(x, p0, obsDays, obsBis, obsOx, EBis, EOx, process, 1);
    problem = createOptimProblem('fmincon', 'x0', x0, 'objective', obj_func, 'Aineq', A, 'bineq', b, 'lb', lb, 'ub', ub, ...
        'options', options);
    % ms = MultiStart('StartPointsToRun', 'bounds', 'TolX', 1e-12, 'TolFun', 1e-12);
    ms = MultiStart('Display', 'off', 'StartPointsToRun', 'all', 'TolX', 1e-12, 'MaxTime', 100);
    
%     gs = GlobalSearch('StartPointsToRun', 'bounds-ineqs', 'TolX', 1e-12);
    
    xminm = [];
     
    while (numel(xminm) == 0)
        [xminm, fmin, ~, output, manyminsm] = run(ms, problem, 10);        
    end
    
    %call DSHydroxy with the optimal parameter values once more to get the hessian 
        %and the model prediction for each time point
    [~, derlogLik, pBis, pOx, pBisModel, pOxModel, pAllStates] = DSHydroxyEmbryoGW(xminm, p0, obsDays, obsBis, obsOx, ...
        EBis, EOx, process, 1);
    H = hessian;
    
    S = diag(diag(H));
    H_scaled = sqrt(S) * H * sqrt(S);
    cond_number = abs(max(eig(H_scaled))/ min(eig(H_scaled)));
    
    %covariance estimate as OPG
%     covEst = inv(-derlogLik' * (-derlogLik));
    
    %if previous covariance is not positive semi-definite
    %-> has negative diagonal elements compute matlab returned hessian
    %close to the argmin point
    if (~isPositiveDefinite(H))
        
%         %fix recognition prob to be 1 (reduce the dimension size by 1 param and re run fmincon)
        lb = [0, -1/lastDay, 0, -1/lastDay, 0, -1/lastDay, 1];
        ub = [1, 1/lastDay, 1, 1/lastDay, 1, 1/lastDay, 1];
        %call fmincon again to compute and return the hessian matrix by matlab
        x0(end) = 1;

        [xminm, ~, ~, ~, ~, ~, ~] = fmincon(obj_func, x0, A, b, [], [], lb, ub, [], options);
        
        %call DSHydroxy with the optimal parameter values once more to get the hessian 
        %and the model prediction for each time point
        [~, ~, pBis, pOx, pBisModel, pOxModel, pAllStates] = DSHydroxyEmbryoGW(xminm, p0, obsDays, obsBis, obsOx, ...
            EBis, EOx, process, 1);
        H = hessian;
        
%         if it is still non positive semi-definite
%         compute the nearest positive semi definite matrix
        index = find(abs(xminm) <= 0.001);
        if (~isPositiveDefinite(H) && ~isempty(index))
            
            if strcmp(bin, '3629073')
                disp('here');
            end    
            %the non spd hessian is attributed to the points which are
            %close to 0. Find them and fix them to 0
            lb(index) = 0;
            ub(index) = 0;
            x0(index) = 0;
            %get new xminm
            [xminm, ~, ~, ~, ~, ~, H] = fmincon(obj_func, x0, A, b, [], [], lb, ub, [], options);
            %get additional indexes smaller than 0.001
            index = setdiff(find(abs(xminm) <= 0.001), index);
            
        end
        
        %if it is still non positive semi-definite
        %compute the nearest positive semi definite matrix
        if (~isPositiveDefinite(H))
            SPDFlag = 1;           
            H = topdm(H);
            
        end

        del = find((ub == lb));
        keep = find((ub ~= lb));
        %throw out of the hessian the nuisance parameters 
        H(del,:) = [];
        H(:,del) = [];
        
             
    end
    
    %square root of the Fisher Information (of the observations) is an
    %approximate lower bound for the standard deviations of the observations
    Cov = inv(H);
        
    %impose missing symmetry from negligible numerical differences
    Cov = makeCovSymmetric(Cov);
    
    
    %throw out of params nan variables
    paramsOut = xminm;
    paramsOut(del) = nan;
    %create CovOut putting nan in entries of fixed parameters
    %and the Cov entry in entries of non-fixed params
    CovOut = nan(numOfParams, numOfParams);
    CovOut(keep, keep) = Cov;
    sigma = (sqrt(diag(CovOut)))';
    
    %volume of the hyperllipse
    confLevel = 0.95;
    volume = normalHyperEllipseVolume(Cov, confLevel);
    
      
    %count all local optimum different than the global one
    localButNotGlobal = 0;
    for i=1:size(manyminsm, 2)
        if(~isequal(round(manyminsm(1,1).X, 2), round(manyminsm(1,i).X, 2)))
            localButNotGlobal = localButNotGlobal + 1;
        end
    end

     
else
    for i=1:numOfDays
        %display and plot the normalized level  taken from the data of all states over time
        %using MLEstimator

        [pAllStates(dayInd(i),:), pBis(dayInd(i),:), pOx(dayInd(i),:), pBisModel(dayInd(i),:), pOxModel(dayInd(i),:)] = ...
            maxLikelihoodGW(obsBis(dayInd(i),:), obsOx(dayInd(i),:), ...
            EBis(:,:,dayInd(i)), EOx(:,:,dayInd(i)), dataType);

    end
            
end

%Kullback Leibler Divergence
KLBis = kullbackLeibler(pBis, pBisModel);
KLOx = kullbackLeibler(pOx, pOxModel);


if ~isempty(chrDir)
    
    %THIS PATH HAS TO BE GIVEN AS AN ARGUMENT BY THE USER!
    %make a directory for that chromosome if it does not exist
    chrAllCpGsDir = strcat(chrDir, '/allCpGs');
    if ~(exist(chrAllCpGsDir, 'dir') == 7)
        mkdir(chrAllCpGsDir);
    end
    outFilePath = strcat(chrAllCpGsDir, '/', stringCpGs, '.txt');
    %prints the filepath in the screen
%     fprintf('%s  ', outFilePath);
    fileID = fopen(outFilePath, 'w', 'ieee-le','UTF-8');

    stringDays = sprintf('%d ', obsDays);

    fprintf(fileID, ['Results for the %dth CpG region of chr%s. It contains the CpGs %s\n with '... 
       'a total number of %d observations at days %s \n'], ...
        bin, chr, stringCpGs, numOfSamples, stringDays);

    fprintf(fileID, 'BS: \t');
    for t=1:numOfDays
        fprintf(fileID, 'd%d %d\t', obsDays(t), numOfSamplesBS(dayInd(t)));
    end
    fprintf(fileID, '\n');
    fprintf(fileID, 'oxBS: \t');
    for t=1:numOfDays
        fprintf(fileID, 'd%d %d\t', obsDays(t), numOfSamplesoxBS(dayInd(t)));
    end
    fprintf(fileID, '\n\n -----------------------------Maximum Likelihood Estimator-------------------------------\n');



    %efficiencies have also been computed from the above if conditions
    if numOfDays == 3
        fprintf(fileID, '%s', output.message);
        fprintf(fileID, ' The local but not global maxima are: %d \n\n', localButNotGlobal); 

        if SPDFlag
            fprintf(fileID, ['Warning: covariances are approximately computed by the inv. of the nearest SPD \n ' ...
                'of the hessian matrix because the numerical hessian was not positive semi-definite \n\n']);
        end

        %Print the estimated parameters the stds
        %and relative half widths
        fprintf(fileID, 'The parameters of the model are: \n\n');
        formatSpec = ['%s = %4.4f ', char(177), ' (%4.4f)\n'];
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

        %print lambda (total methylation) coef
        lambdaCoef(paramsOut, CovOut, fileID);

        %print the hypervolume of the confidence region 
        fprintf(fileID, 'The volume of the hyperellipse of the multivariate normal is %d (confLevel = %.2f)\n\n', ...
            volume, confLevel);



    end

    fprintf(fileID, '---------------------------------------------------------------\n');

    C = num2cell(pBis(dayInd, :));
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

    fprintf(fileID, 'The data distribution of BS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec, C{row,:});
    end
    fprintf(fileID, '\n');

    C = num2cell(pBisModel(dayInd,:));
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

    fprintf(fileID, 'The predicted distribution of BS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec, C{row,:});
    end
    fprintf(fileID, '\n');
    fprintf(fileID, 'The KL divergence for BS prediction is %.4f\n\n', KLBis);

    fprintf(fileID, '---------------------------------------------------------------\n');

    C = num2cell(pOx(dayInd,:));
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

    fprintf(fileID, 'The data distribution of oxBS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec, C{row,:});
    end
    fprintf(fileID, '\n');

    C = num2cell(pOxModel(dayInd,:));
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

    fprintf(fileID, 'The predicted distribution of oxBS states \n\n');
    fprintf(fileID, 'data TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec, C{row,:});
    end
    fprintf(fileID, '\n');
    fprintf(fileID, 'The KL divergence for oxBS prediction is %.4f \n\n', KLOx);

    fprintf(fileID, '---------------------------------------------------------------\n');
    KL = KLBis + KLOx;
    fprintf(fileID, 'The total KL divergence for MLE is %1.4f\n\n', KL);
    fprintf(fileID, '---------------------------------------------------------------\n');

    C = num2cell(pAllStates(dayInd,:));
    C = [dataLabels, C];
    formatSpec = '%s %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n';

    fprintf(fileID, 'The predicted distribution of the hidden states \n\n');
    fprintf(fileID, 'data uu um mu uh hu hm mh mm hh \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, formatSpec, C{row,:});
    end

    fclose(fileID);
 
end

end
