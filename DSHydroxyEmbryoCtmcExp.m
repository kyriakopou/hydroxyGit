function [nL, der, pData, pModel, pAllStates] = DSHydroxyEmbryoCtmcExp(x, maxDegree, p0, days, obsMatrix, EMatrix, ~, derCompFlag)
%EIAP_DSHYDROXY summary of this function goes here
%   Detailed explanation goes here

% x(1): maintenance b0
% x(2): maintenance b1
% x(3): de-Novo b0
% x(4): de-Novo b1
% x(5): hydroxylation b0
% x(6): hydroxylation b1 
% x(7): demethylation b0
% x(8): demethylation b1
% x(9): prob h is recognized as u

global Hessian

if (numel(x) == 0)
    disp('numel(x) = 0');
    return;
% elseif any(A * x' > b);
%     disp('negative effic. over time');
% %     return;
end


% constants
numOfDataPoints = size(days, 1);
numObsStates = size(obsMatrix, 2);
numOfParams = size(x, 2);
numOfExperiments = size(obsMatrix, 3);
numOfStates = getNumberOfStatesTemp(numOfExperiments);

%modify process vector to get entries 
%for all dataPoints
% processAll = cell(days(end), 1);
% 
% for i=2:numOfDataPoints
%     processAll{days(i)} = process(i);
% end
% for i=days(end)-1:-1:1
%     if (isempty(processAll{i}))
%         processAll{i} = processAll{i+1};
%     end
% end    



%normalized Observations for all time points and all experiments
pData = obsMatrix ./ repmat(sum(obsMatrix, 2), 1, numObsStates);


%------Estimate the transient probabilites for each observable day
%Bisulfite - OxBisulfite
%the estimation for p0 are the normalized frequencies of the observations
pModel = zeros(numOfDataPoints, numObsStates, numOfExperiments);
for e=1:numOfExperiments
    pModel(1,:,e) = pData(1,:,e);
end

pAllStates = zeros(numOfDataPoints, numOfStates);
%the derivatives matrices
derpModel = zeros(numOfDataPoints, numObsStates, numOfParams, numOfExperiments);
%second partial derivatives matrices
secDerpModel = zeros(numOfDataPoints, numObsStates, numOfParams, numOfParams, numOfExperiments);
%the derivatives temporary vectors derpPrev, derpNext
derpAllStates = zeros(numOfDataPoints, numOfStates, numOfParams);
%the second partial derivatives temporary vectors 
secDerpAllStates = zeros(numOfDataPoints, numOfStates, numOfParams, numOfParams);


%compute the probability mass pNext and the derivative prob mass 
%for the next day. In days of observations weightening both with Cbis, COxbis
k = 2;
pAllStates(1,:) = p0;
p = p0;
% Create constant processes matrices


for t=1:days(end)
    %mu2 and matrices M, derM change at each day bec. mu2 changes
    [mu1, mu2, h, f, dem] = getEfficiencies(x, maxDegree, t);
    p_rec = x(end);
    
    
    %update the p vector;
    
    %get the Q matrix and initialize its derivative to 0s (t should be included to the arg list for the derivative)
    [Q, derQ, secondDerQ] = QMatrixAndDerivatives(dem, f, h, mu2, t);
    derExpQd = zeros(size(derQ));
    
    p = p * myCellDivision(0.5) * myMaintenanceAct(mu1, p_rec, t) * expm(Q*t);
    
    %if d = observation day compute the prob. of each of the observable
    %states for Bis and OxBis
    if t == days(k)
        
        pAllStates(k,:) = p;
        
        for e=1:numOfExperiments
            pModel(k,:,e) = pAllStates(k,:) * EMatrix(:,:,k,e);
        end
        
        if derCompFlag
            for i=1:numOfParams
                [~, derExpQd(:,:,i)] = dexpmt(Q*t, t*derQ(:,:,i));
                derpAllStates(k,:,i) = p0*derExpQd(:,:,i);
                
                for e=1:numOfExperiments    
                    derpModel(k,:,i,e) = derpAllStates(k,:,i) * EMatrix(:,:,k,e);
                end
                for j=1:numOfParams
%                     secDerExpQd = t*(secondDerQ(:,:,i,j)*expm(Q*t) + derQ(:,:,j)*derExpQd(:,:,i));
                    secDerExpQd = zeros(size(Q));
                    secDerpAllStates(k,:,i,j) = p0*secDerExpQd; 
                    for e=1:numOfExperiments
                        secDerpModel(k,:,i,j,e) = secDerpAllStates(k,:,i,j) * EMatrix(:,:,k,e);
                    end
                end    
            end  
        end    
        k = k+1;     
    end
    
end


%compute likelihood summing over all states over all time points and over
%all experiments
logLik = nansum(nansum(nansum(obsMatrix .* log(pModel))));

%First-second derivatives of the log-likelihood
derlogLik = zeros(1, numOfParams);
secDerlogLik = zeros(numOfParams, numOfParams);
if derCompFlag
    for i=1:numOfParams
        %reshape dimensions of der matrix
        derpMTempI = reshape(derpModel(:,:,i,:), [size(derpModel(:,:,i,:), 1), size(derpModel(:,:,i,:), 2), size(derpModel(:,:,i,:), 4)]); 
        derlogLik(i) = nansum(nansum(nansum(obsMatrix .* derpMTempI ./ pModel)));
        %second derivatives of the log-likelihood
        for j=1:numOfParams
            %reshape dimensions of der and secDer matrices
            derpMTempJ = reshape(derpModel(:,:,j,:), [size(derpModel(:,:,j,:), 1), size(derpModel(:,:,j,:), 2), size(derpModel(:,:,j,:), 4)]); 
            secDerpMTempIJ = reshape(secDerpModel(:,:,i,j,:), [size(secDerpModel(:,:,i,j,:), 1), size(secDerpModel(:,:,i,j,:), 2), size(secDerpModel(:,:,i,j,:), 5)]); 
            secDerlogLik(i,j) = nansum(nansum(nansum(obsMatrix .* (secDerpMTempIJ .* pModel - derpMTempI .* derpMTempJ) ./ pModel.^2 )));
        end    
    end
end

%transform to a minimization problem
nL = -logLik;
der = -derlogLik;
der = [0 0 der 0];
%update the hessian in the global variable
Hessian = -secDerlogLik;


end

