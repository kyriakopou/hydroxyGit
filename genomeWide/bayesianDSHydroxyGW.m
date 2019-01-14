%parameter estimation procedure
%CALL THIS FUNCTION INSTEAD OF GUI WHEN YOU WANT RESULTS FOR MANY REGIONS
%(either singleCpGs or whole regions) -- TO STORE INTO FIXED FOLDER etc.
function [meanPosterior, CovPosterior, pAllStates, KLBis, KLOx, volume, accRate] = bayesianDSHydroxyGW(binAndCpGs, ...
    maxObsDaysVector, maxDataLabels, obsBis, obsOx, errorsBis, errorsOx, process, p0, meanPrior, CovPrior, ...
    paramNames, chrDir)

%days with non zeros observations
ind = find(~all(obsBis == 0, 2));
obsDays = maxObsDaysVector(ind);
dataLabels = maxDataLabels(ind);

numOfDays = size(obsDays, 1);
numObsStates = size(obsBis, 2);
numOfStates = 9;
numOfParams = size(meanPrior, 2);
CovPosterior = nan(numOfParams, numOfParams);


bin = cell2mat(binAndCpGs(1));
cpgArray = cell2mat(binAndCpGs(2:end));
stringCpGs = sprintf('%d_', cpgArray);
regionName = strcat('bin_', int2str(bin));
dataType = 'oxBS';


% CONVERSION ERRORS
%Conversion Matrix bisulfite - Oxbisulfite for each day
EBis = zeros(numOfStates, numObsStates, numOfDays);
EOx = zeros(numOfStates, numObsStates, numOfDays);
for t=1:numOfDays
    EBis(:,:,t) = convErrorBis(errorsBis(t,:));
    EOx(:,:,t) = convErrorOx(errorsOx(t,:));   
end

lastDay = obsDays(numOfDays);

CovPriorDef = ...
[   0.3947   -0.0679   -0.1393    0.0273    0.1416   -0.0339    0.0601
   -0.0679    0.0162    0.0225   -0.0041   -0.0227    0.0070    0.0001
   -0.1393    0.0225    0.0732   -0.0141   -0.0491    0.0106   -0.0007
    0.0273   -0.0041   -0.0141    0.0029    0.0091   -0.0019    0.0035
    0.1416   -0.0227   -0.0491    0.0091    0.0815   -0.0193   -0.0123
   -0.0339    0.0070    0.0106   -0.0019   -0.0193    0.0058    0.0019
    0.0601    0.0001   -0.0007    0.0035   -0.0123    0.0019    0.1657];

%edit the prob entry if it is close to 0 or nan
nanParamEntries = find(isnan(meanPrior));
meanPrior(nanParamEntries) = 0;
%if rec prob is nan or very small
if ismember(7, nanParamEntries) || meanPrior(end) < 0.1
    meanPrior(end) = 0.75;
end

nanParamOrLargeStdEntries = find(isnan(CovPrior) | abs(CovPrior) > 100);
CovPrior(nanParamOrLargeStdEntries) = CovPriorDef(nanParamOrLargeStdEntries);

x0 = meanPrior;
sigmas = sqrt(diag(CovPrior))';

%confidence Intervals for meanPrior
confIntLB = meanPrior - 1.96*sigmas;
confIntUB = meanPrior + 1.96*sigmas;

%define the parameter bounds (Initialize both to the same value to fix a parameter)
%linear maintenance, de novo, hydroxy
lb = [repmat([0, -1/lastDay], [1, (numOfParams-1)/2]), 0];
ub = [repmat([1, 1/lastDay], [1, (numOfParams-1)/2]), 1];

%get const entries of feasible region
lbConst = lb(1:2:numOfParams);
ubConst = ub(1:2:numOfParams);


%compute the intersection between the confidence
%intervals of mean prior and the feasible region
confIntLB = max(confIntLB, lb);
confIntUB = min(confIntUB, ub);


% [x0, ~, ~] = unifDependRnd(lbConst, ubConst, lastDay);


%-----METROPOLIS HASTINGS-----
numOfUpdates = 1;
nsamples = 10^4;
confLevel = 0.95;
jointPosterSamples = zeros(nsamples, numOfParams, numOfUpdates);
% priorPdf = @(x) prod(unifpdf(x, confIntLB, confIntUB));
% priorPdf = @(x) unifDependPdf(x, lastDay);
priorPdf = @(x) myTruncNormalPdf(x, meanPrior, sigmas, lastDay);
%likelihood pdf
loglikelihoodPdf = @(x) -DSHydroxyEmbryoGW(x, p0, obsDays, obsBis, obsOx, EBis, EOx, process, 0);

%target distribution
logPdf = @(x) computeTargetPdf(x, priorPdf, loglikelihoodPdf);
    % proposal density and proposal random number generator
%     proppdf = @(x,y) mvnpdf(y, x, ones(1, numOfParams/75));
%     proppdf = @(x,y) unifpdf(y, finalIntLB, finalIntUB);
%     proprnd = @(x) unifrnd(finalIntLB, finalIntUB);

%scale sigma for "optimal" proposal distribution variance
sigmasProp =  sigmas / 7.5;

for k=1:numOfUpdates
    rng default  % For reproducibility
    
    accRate = 0;
    while accRate == 0
     
        proppdf = @(y,x) myTruncNormalPdf(y, x, sigmasProp, lastDay);
        %WE HAVE TO GET A PROPOSAL DIS WITH SMALLER VARIANCE (TO INCREASE THE ACC RATE)

    %     proprnd = @(x) unifDependRnd(lbConst, ubConst, lastDay);
    %     proprnd = @(x) mvnrnd(meanPrior, CovPrior);
    %     proprnd = @(x) mvnrnd(x, ones(1, numOfParams)/75);

        proprnd = @(x) myTruncNormalRnd(x, sigmasProp, lastDay);            

        [jointPosterSamples(:,:,k), accRate] = mhsample(x0, nsamples, 'logpdf', logPdf, 'proprnd', proprnd, ...
            'proppdf', proppdf, 'burnin', 500);
    %     [jointPosterSamples(:,:,k), accRate] = mhsample(x0, nsamples, 'logpdf', logPdf, 'symmetric', true, ...
    %'proprnd', proprnd);

        sigmasProp = sigmasProp / 2;
    
    end
    
    %mean and var of marginal posteriors
    meanPosterior = mean(jointPosterSamples(:,:,k));
    if accRate ~= 0
        CovPosterior = cov(jointPosterSamples(:,:,k));
    end
    stdsPosterior = sqrt(diag(CovPosterior))';
    
    %plot the sampled marginal posteriors for the params
%     fig = figure('Visible', 'off');
    credIntervalsFlag = 0;
    if credIntervalsFlag 
        for i=1:numOfParams
    %         subplot(4,2,i)

    %         h = histfit(jointPosterSamples(:,i,k), 10, 'kernel');
            %normalize the histogram values (both fitted pdf and hist)
    %         h(1).FaceColor = [.8 .8 1];
    %         h(1).YData = h(1).YData / nsamples;
    %         h(2).YData = h(2).YData / nsamples;


            [f, xi, bw] = ksdensity(jointPosterSamples(:,i,k), 'Support', [lb(i), ub(i)]);
            %smoothen the fitted kernel by doubling the bandwidth
            [f, xi] = ksdensity(jointPosterSamples(:,i,k), 'bandwidth', 3*bw, 'Support', [lb(i), ub(i)]);

    %         hold on;
    %         plot(xi, f);

            %CHECKING: fitted dist mean and var
    %         trapz(xi, f)
    %         trapz(xi, f.*xi)

            %find the closest elm of xi to meanPrior(i) 
            [~, elmMLE] = min(abs(xi-meanPrior(i)));

            %plot mle estimator and confidence intervals
    %         line([meanPrior(i), meanPrior(i)], get(gca, 'ylim'), 'Color', 'red', 'LineStyle', ':', 'LineWidth', 1);
    %         line([confIntLB(i), confIntLB(i)], get(gca, 'ylim')/10, 'Color', 'red', 'LineWidth', 1);
    %         line([confIntUB(i), confIntUB(i)], get(gca, 'ylim')/10, 'Color', 'red', 'LineWidth', 1);
    %         
            %find the closest elm of xi to meanPosterior(i) 
            [~, elmBI] = min(abs(xi-meanPosterior(i)));

            %plot BI posterior mean and credible intervals
    %         line([meanPosterior(i), meanPosterior(i)], get(gca, 'ylim'), 'LineStyle', ':', 'LineWidth', 1);

            %compute HPD (highest posterior density credible set) 
            [fSorted, sortInd] = sort(f, 'descend');
            xiToSortedF = xi(sortInd);
            totSum = sum(fSorted);
            HPDInd = ~(cumsum(fSorted) / totSum > confLevel);
            xiHPD = xiToSortedF(HPDInd);
            xiHPDSorted = sort(xiHPD);

            if size(xiHPDSorted,2) > 1
                credIntLB = max(xiHPDSorted(1), lb(i));
                credIntUB = min(xiHPDSorted(end), ub(i));

    %             line([credIntLB, credIntLB], get(gca, 'ylim')/10, 'LineWidth', 1);
    %             line([credIntUB, credIntUB], get(gca, 'ylim')/10, 'LineWidth', 1);
            end

    %         xlabel(paramNames(i), 'FontWeight', 'bold');
    %         if i==1
    %             l = legend({'posterior', 'MLE'}, 'Location', 'northwest', 'Box', 'off');
    %             set(l, 'FontSize', 6);
    %         end    
        end
    end    
%     REMOVE COMMENTS TO GET PNG FILES
%     %check if the folder for png files is there otherwise create it
%     picPath = strcat(resultsPath, '/', chr, '/', 'posterFigs');
%     if ~(exist(picPath, 'file') == 7)
%         mkdir(picPath);
%     end
%     picPath = strcat(picPath, '/CpG_', stringCpGs);
%     print(fig, picPath, '-dpng');
    
    %update the priorPdf (if more than 1 bayesian repetitions)
    priorPdf = @(x) mvnpdf(x, meanPosterior, CovPosterior);
end    

%get the output of the model for the means of the sampled parameters' posteriors
[~, ~, pBis, pOx, pBisModel, pOxModel, pAllStates] = DSHydroxyEmbryoGW(meanPosterior, p0, obsDays, obsBis, obsOx, ...
    EBis, EOx, process, 0);

%get the volume of the hyperellipse of the normal dist
volume = normalHyperEllipseVolume(CovPosterior, confLevel);

%Kullback Leibler Divergence 
KLBis = kullbackLeibler(pBis, pBisModel);
KLOx = kullbackLeibler(pOx, pOxModel);

if ~isempty(chrDir)
    
    %THIS PATH HAS TO BE GIVEN AS AN ARGUMENT BY THE USER!
    chrAllCpGsDir = strcat(chrDir, '/allCpGs');
    outFilePath = strcat(chrAllCpGsDir, '/', stringCpGs, '.txt');

    fileID = fopen(outFilePath, 'a', 'ieee-le','UTF-8');

    %efficiencies have also been computed from the above if conditions
    if numOfDays == 3

        %replace with bayesian values
        paramsOut = meanPosterior;
        Cov = CovPosterior;
        sigmas = stdsPosterior;

        fprintf(fileID, '\n\n -----------------------------Bayesian Inference------------------------------\n\n');
        %Print the estimated parameters the stds
        %and relative half widths
        fprintf(fileID, 'The parameters of the model (acc. rate was %.2f) are: \n\n', accRate);
        formatSpec = ['%s = %4.4f ', char(177), ' (%4.4f)\n'];
        for i = 1:size(paramsOut, 2)
            %fprintf(fileID, formatSpec, paramNames{1,:});
            %formatSpec = '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f\n';
            fprintf(fileID, formatSpec, char(paramNames(i)), paramsOut(i), sigmas(i));
        end    
        fprintf(fileID, '\n');

        lambdaCoef(meanPosterior, Cov, fileID)
        %print the hypervolume of the confidence region 
        fprintf(fileID, 'The volume of the hyperellipse of the multivariate normal is %d (confLevel = %.2f)\n', ...
            volume, confLevel);

    end

        fprintf(fileID, '---------------------------------------------------------------\n');


        C = num2cell(pBis);
        C = [dataLabels, C];
        formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The data distribution of BS states \n\n');
        fprintf(fileID, 'data TT TC CT CC \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, formatSpec, C{row,:});
        end
        fprintf(fileID, '\n');

        C = num2cell(pBisModel);
        C = [dataLabels, C];
        formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The predicted distribution of BS states \n\n');
        fprintf(fileID, 'data TT TC CT CC \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, formatSpec, C{row,:});
        end
        fprintf(fileID, '\n');
        fprintf(fileID, 'The KL divergence for BS prediction is %.4f \n\n', KLBis);

        fprintf(fileID, '---------------------------------------------------------------\n');

        C = num2cell(pOx);
        C = [dataLabels, C];
        formatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The data distribution of oxBS states \n\n');
        fprintf(fileID, 'data TT TC CT CC \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, formatSpec, C{row,:});
        end
        fprintf(fileID, '\n');

        C = num2cell(pOxModel);
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
        fprintf(fileID, 'The total KL divergence for BI is %1.4f \n\n', KL);
        fprintf(fileID, '---------------------------------------------------------------\n');

        C = num2cell(pAllStates);
        C = [dataLabels, C];
        formatSpec = '%s %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The predicted distribution of the hidden states \n\n');
        fprintf(fileID, 'data uu um mu uh hu hm mh mm hh \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, formatSpec, C{row,:});
        end
        fprintf(fileID, '\n');

    fclose(fileID);    
end    

end
