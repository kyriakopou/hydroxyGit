function [nL, der, pData, pModel, pAllStates] = DSHydroxyEmbryo_ctmc(x, maxDegree, knots, p0, days, ...
    obsMatrix, EMatrix, process, derCompFlag, matrixFileName)
%EIAP_DSHYDROXY summary of this function goes here
%   Detailed explanation goes here


global Hessian

if (numel(x) == 0)
    disp('numberOfElements of x = 0');
    return;
% elseif any(A * x' > b);
%     disp('negative effic. over time');
% %     return;
end


% constants
numOfDataPoints = length(days);
numObsStates = size(obsMatrix, 2);
numOfParams = size(x, 2);
numOfExperiments = size(obsMatrix, 3);

numOfStates = getNumberOfStatesTemp(numOfExperiments);

%to get entries for all dataPoints
processAll = getProcessForEachDay(days, process);

%normalized Observations for all time points and all experiments
pData = obsMatrix ./ repmat(sum(obsMatrix, 2), 1, numObsStates);

%------Estimate the transient probabilites for each observable day
%Bisulfite - OxBisulfite
%the estimation for pModel at time point 0
pModel = zeros(numOfDataPoints, numObsStates, numOfExperiments);
for e=1:numOfExperiments
    pModel(1,:,e) = p0 * EMatrix(:,:,1,e);
end

pAllStates = [p0; zeros(numOfDataPoints-1, numOfStates)];
%the derivatives matrices
derpModel = zeros(numOfDataPoints, numObsStates, numOfParams, numOfExperiments);
%second partial derivatives matrices
secDerpModel = zeros(numOfDataPoints, numObsStates, numOfParams, numOfParams, numOfExperiments);
%the derivatives temporary vectors derpPrev, derpNext
derpAllStates = zeros(numOfDataPoints, numOfStates, numOfParams);



%compute the probability mass pNext and the derivative prob mass 
%for the next day. In days of observations weightening both with Cbis, COxbis
% Create constant processes matrices
k=2;

%initialize prob and der vector
switch derCompFlag
    case 2
        pAndDerPrev = [p0, zeros(1, numOfParams * numOfStates), zeros(1, numOfParams^2 * numOfStates)];
    case 1
        pAndDerPrev = [p0, zeros(1, numOfParams * numOfStates)];
    case 0
        pAndDerPrev = p0;
end   

for t=1:days(end)
    
    EFF = getEfficiencies(x, maxDegree, t, knots);
    
    if (strcmp(processAll{t}, 'rep'))  
        %------this block executes only if we have cell division
        mu1 = EFF(:,1);
        p_rec = x(end);
        maintFunction = str2func(matrixFileName{1});
        
%         [M_maint, derM_maint, secDerM_maint] = myMaintenanceAct(mu1, p_rec, t);
        [M_maint, derM_maint, secDerM_maint] = maintFunction(mu1, p_rec, t);

        pVector = pAndDerPrev(1:numOfStates);
        firstDer = pAndDerPrev(numOfStates+1:(numOfParams+1)*numOfStates);
        secDer = pAndDerPrev((numOfParams+1)*numOfStates:end);

        CD = myCellDivision(0.5);
        %if derivative computation is on update derivative vectors
        if derCompFlag >=1
            
            %the derivatives temporary vectors derpPrev, derpNext
            derp = zeros(numOfParams, numOfStates);
            secDerp = zeros(1, numOfStates, numOfParams, numOfParams);
          
            for i=1:numOfParams
                derp(i,:) = firstDer((i-1)*numOfStates+1:i*numOfStates) * CD * M_maint + ...
                    pVector * CD * derM_maint(:,:,i);
                if derCompFlag == 2
                    %update derivative parts
                    for j=1:numOfParams
                        secDerp(:,:,i,j) = secDer((i-1)*numOfParams+(j-1)*numOfStates+1:(i-1)*numOfParams+j*numOfStates) ...
                            * CD * M_maint + firstDer((j-1)*numOfStates+1:j*numOfStates) * CD * derM_maint(:,:,i) + ...
                            firstDer((i-1)*numOfStates+1:i*numOfStates) * CD * derM_maint(:,:,j) + ...
                            pVector * CD * secDerM_maint(:,:,i,j);
                    end    
                    %secDer
                    pAndDerPrev((numOfParams+1)*numOfStates+1:end) = reshape(secDerp, [numOfParams^2 * numOfStates, 1]);
                end    
                %firstDer
                pAndDerPrev(numOfStates+1:(numOfParams+1)*numOfStates) = reshape(derp', [numOfStates * numOfParams, 1]);

            end                      
        end
        %update probability vector
        pAndDerPrev(1:numOfStates) = pAndDerPrev(1:numOfStates) * CD * M_maint;
    end
    
    %------------
%     call the ode-solver to get the solution within        
%     options = odeset('RelTol', 1e-2, 'AbsTol', 1e-5); 
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6); %default values rel = 1e-3, abs = 1e-6
    QMatrixFunction = str2func(matrixFileName{end});
    [T, pAndDerNext] = ode45(@(time, p) masterEquation(time, p, x, maxDegree, knots, numOfStates, derCompFlag, ...
        QMatrixFunction), [t-1, t], pAndDerPrev, options);
    
    %get next transient and derivative probability vector pAndDer
    pAndDerPrev = pAndDerNext(end,:);
      
    %store it to pAllStates if t is an observation time day
    if t == days(k)
%         ind = find(T > t, 1, 'first');
%         pAndDerNext(end,:)
        
        pAllStates(k,:) = pAndDerPrev(1:numOfStates);
        
        %pModel is the dist of the obs. States for each experiment
        for e=1:numOfExperiments
            pModel(k,:,e) = pAllStates(k,:) * EMatrix(:,:,k,e);
        end
                
        
        if derCompFlag
            for i=1:numOfParams
                
                derpAllStates(k,:,i) = pAndDerPrev(numOfStates*i+1:numOfStates*(i+1));

                for e=1:numOfExperiments    
                    derpModel(k,:,i,e) = derpAllStates(k,:,i) * EMatrix(:,:,k,e);
                end
            end
            
            if derCompFlag == 2
                secDer = pAndDerPrev((numOfParams+1)*numOfStates+1:end);
                for j=1:numOfParams
                    for e=1:numOfExperiments
                        IJ_indStart = (i-1)*numOfParams+(j-1)*numOfStates+1;
                        IJ_indEnd = (i-1)*numOfParams+j*numOfStates;
                        secDerpModel(k,:,i,j,e) = secDer(IJ_indStart: IJ_indEnd) * EMatrix(:,:,k,e);
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
if derCompFlag >=1
    for i=1:numOfParams
        %reshape dimensions of der matrix
        derpMTempI = reshape(derpModel(:,:,i,:), [size(derpModel(:,:,i,:), 1), size(derpModel(:,:,i,:), 2), ...
            size(derpModel(:,:,i,:), 4)]); 
        derlogLik(i) = nansum(nansum(nansum(obsMatrix .* derpMTempI ./ pModel)));
        if derCompFlag == 2
            %second derivatives of the log-likelihood
            for j=1:numOfParams
                %reshape dimensions of der and secDer matrices
                derpMTempJ = reshape(derpModel(:,:,j,:), [size(derpModel(:,:,j,:), 1), size(derpModel(:,:,j,:), 2), ...
                    size(derpModel(:,:,j,:), 4)]); 
                secDerpMTempIJ = reshape(secDerpModel(:,:,i,j,:), [size(secDerpModel(:,:,i,j,:), 1), ...
                    size(secDerpModel(:,:,i,j,:), 2), size(secDerpModel(:,:,i,j,:), 5)]); 
                secDerlogLik(i,j) = nansum(nansum(nansum(obsMatrix .* (secDerpMTempIJ.*pModel - derpMTempI.*derpMTempJ) ...
                    ./ pModel.^2 )));
            end   
        end    
    end
end

%transform to a minimization problem
nL = -logLik;
der = -derlogLik;
%update the hessian in the global variable
Hessian = -secDerlogLik;


end

