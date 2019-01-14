function [nL, der, pData, pModel, pAllStates] = DSHydroxyEmbryo_dtmc(x, maxDegree, knots, p0, days, ...
    obsMatrix, EMatrix, process, derCompFlag, matrixFileName)
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
s = 0.5;        % cell division (0.5 represents equally likely strands)

p_u = x(end);   % probability that h is seen as u

numOfTimePoints = size(days, 1);
numObsStates = size(obsMatrix, 2);
numOfParams = size(x, 2);
numOfV = 11;
numOfExperiments = size(obsMatrix, 3);
numOfStates = getNumberOfStatesTemp(numOfExperiments);


processAll = getProcessForEachDay(days, process);



%normalized Observations for all time points and all experiments
pData = obsMatrix ./ repmat(sum(obsMatrix, 2), 1, numObsStates);


%------Estimate the transient probabilites for each observable day
%Bisulfite - OxBisulfite
%the estimation for of pModel at time point 0
pModel = zeros(numOfTimePoints, numObsStates, numOfExperiments);
for e=1:numOfExperiments
    pModel(1,:,e) = p0 * EMatrix(:,:,1,e);
end
pAllStates = zeros(numOfTimePoints, numOfStates);
%the derivatives matrices
derpModel = zeros(numOfTimePoints, numObsStates, numOfParams, numOfExperiments);
%second partial derivatives matrices
secDerpModel = zeros(numOfTimePoints, numObsStates, numOfParams, numOfParams, numOfExperiments);
%the derivatives temporary vectors derpPrev, derpNext
derpPrev = zeros(1, numOfStates, numOfParams);
derpNext = zeros(1, numOfStates, numOfParams);
%derivatives in case p0 is also param
% for i=numOfV+1 : numOfParams
%    derpPrev(:, i-numOfV, i) = 1; 
% end    
%the second partial derivatives temporary vectors 
secondDerpPrev = zeros(1, numOfStates, numOfParams, numOfParams);
secondDerpNext = zeros(1, numOfStates, numOfParams, numOfParams);
%initialize derivatives of P matrix 
derP = zeros(numOfStates, numOfStates, numOfParams);
%initialize second partial derivatives P matrix 
secondDerP = zeros(numOfStates, numOfStates, numOfParams, numOfParams);


%compute the probability mass pNext and the derivative prob mass 
%for the next day. In days of observations weightening both with Cbis, COxbis
k = 2;
pPrev = p0;
pAllStates(1,:) = p0;
% Create constant processes matrices
D = myCellDivision(s);

maintenance = str2func(matrixFileName{1});
deNovo = str2func(matrixFileName{2});
hydroxylation = str2func(matrixFileName{3});
formal = str2func(matrixFileName{4});
active = str2func(matrixFileName{5});

for d=1:days(end)
    
    eff = getEfficiencies(x, maxDegree, d);
    mu1 = eff(:,1);
    mu2 = eff(:,2);
    h = eff(:,3);
    f = eff(:,4);
    dem = eff(:,5);
    
%     [MOld, derMOld, secondDerMOld] = methylation(mu1, mu2, p_u, d);
    
    [M, derM, secondDerM] = maintenance(mu1, p_u, d);
    [Md, derMd, secondDerMd] = deNovo(mu2, d);
    [H, derH, secondDerH] = hydroxylation(h, d);
    [F, derF, secondDerF] = formal(f, d);
    [A, derA, secondDerA] = active(dem, d);
    
    
    %Compute P Matrix for day d
    if (strcmp(processAll{d}, 'rep'))       
        P = D * M * Md * H * F * A;
    elseif (strcmp(processAll{d}, 'no-rep'))
        P = Md * H * F * A; 
    end    
    %compute the derivatives of P matrix wrt to every parameter x_i
    if derCompFlag
        for i=1:numOfV
            if (strcmp(processAll{d}, 'rep'))
                derP(:,:,i) = D * (derM(:,:,i)*Md*H*F*A + M*derMd(:,:,i)*H*F*A + M*Md*derH(:,:,i)*F*A + M*Md*H*derF(:,:,i)*A + M*Md*H*F*derA(:,:,i));
                %compute the second partial derivative (i,j) of P matrix 
                for j=1:numOfV
                    secondDerP(:,:,i,j) = D * ( ...
                    derM(:,:,i)*derMd(:,:,j)*H*F*A + derM(:,:,j)*derMd(:,:,i)*H*F*A + ...
                    derM(:,:,i)*Md*derH(:,:,j)*F*A + derM(:,:,j)*Md*derH(:,:,i)*F*A + ...
                    derM(:,:,i)*Md*H*derF(:,:,j)*A + derM(:,:,j)*Md*H*derF(:,:,i)*A + ...
                    derM(:,:,i)*Md*H*F*derA(:,:,j) + derM(:,:,j)*Md*H*F*derA(:,:,i) + ...
                    M*derMd(:,:,i)*derH(:,:,j)*F*A + M*derMd(:,:,j)*derH(:,:,i)*F*A + ...
                    M*derMd(:,:,i)*H*derF(:,:,j)*A + M*derMd(:,:,j)*H*derF(:,:,i)*A + ...
                    M*derMd(:,:,i)*H*F*derA(:,:,j) + M*derMd(:,:,j)*H*F*derA(:,:,i) + ... 
                    M*Md*derH(:,:,i)*derF(:,:,j)*A + M*Md*derH(:,:,j)*derF(:,:,i)*A + ...
                    M*Md*derH(:,:,i)*F*derA(:,:,j) + M*Md*derH(:,:,j)*F*derA(:,:,i) + ...
                    M*Md*H*derF(:,:,i)*derA(:,:,j) + M*Md*H*derF(:,:,j)*derA(:,:,i) + ...
                    secondDerM(:,:,i,j)*Md*H*F*A + M*secondDerMd(:,:,i,j)*H*F*A + M*Md*secondDerH(:,:,i,j)*F*A + ...
                    M*Md*H*secondDerF(:,:,i,j)*A + M*Md*H*F*secondDerA(:,:,i,j));    
                
%                     secondDerP(:,:,i,j) = D * ( ...
%                     derM(:,:,i)*derMd(:,:,j)*H*F + derM(:,:,j)*derMd(:,:,i)*H*F + ...
%                     derM(:,:,i)*Md*derH(:,:,j)*F + derM(:,:,j)*Md*derH(:,:,i)*F + ...
%                     derM(:,:,i)*Md*H*derF(:,:,j) + derM(:,:,j)*Md*H*derF(:,:,i) + ...
%                     M*derMd(:,:,i)*derH(:,:,j)*F + M*derMd(:,:,j)*derH(:,:,i)*F + ...
%                     M*derMd(:,:,i)*H*derF(:,:,j) + M*derMd(:,:,j)*H*derF(:,:,i) + ...
%                     M*Md*derH(:,:,i)*derF(:,:,j) + M*Md*derH(:,:,j)*derF(:,:,i) + ...
%                     secondDerM(:,:,i,j)*Md*H*F + M*secondDerMd(:,:,i,j)*H*F + M*Md*secondDerH(:,:,i,j)*F + M*Md*H*secondDerF(:,:,i,j));   
                
                end                               
            elseif (strcmp(processAll{d}, 'no-rep'))
                derP(:,:,i) = derMd(:,:,i)*H*F*A + Md*derH(:,:,i)*F*A + Md*H*derF(:,:,i)*A + Md*H*F*derA(:,:,i);
                for j=1:numOfV
                    secondDerP(:,:,i,j) = ...
                    derMd(:,:,i)*derH(:,:,j)*F*A + derMd(:,:,j)*derH(:,:,i)*F*A +...
                    derMd(:,:,i)*H*derF(:,:,j)*A + derMd(:,:,j)*H*derF(:,:,i)*A + ...
                    derMd(:,:,i)*H*F*derA(:,:,j) + derMd(:,:,j)*H*F*derA(:,:,i) + ...
                    Md*derH(:,:,i)*derF(:,:,j)*A + Md*derH(:,:,j)*derF(:,:,i)*A + ...
                    Md*derH(:,:,i)*F*derA(:,:,j) + Md*derH(:,:,j)*F*derA(:,:,i) + ...
                    Md*H*derF(:,:,i)*derA(:,:,j) + Md*H*derF(:,:,j)*derA(:,:,i) + ...
                    secondDerMd(:,:,i,j)*H*F*A + Md*secondDerH(:,:,i,j)*F*A + Md*H*secondDerF(:,:,i,j)*A + Md*H*F*secondDerA(:,:,i,j);
                end
            end
        end
    end
    %compute the prob. mass at day d
    pNext = pPrev * P;  
    
    if derCompFlag 
        %compute the derivatives wrt to all parameters  
        for i=1:numOfParams
                derpNext(:,:,i) = derpPrev(:,:,i) * P + pPrev * derP(:,:,i);
                %compute the second derivatives
                for j=1:numOfParams
                    secondDerpNext(:,:,i,j) = secondDerpPrev(:,:,i,j) * P + derpPrev(:,:,i) * derP(:,:,j) + ...
                    derpPrev(:,:,j) * derP(:,:,i) + pPrev * secondDerP(:,:,i,j) ;
                end 
        end
    end
   
    %if d = observation day compute the prob. of each of the observable
    %states for Bis and OxBis
    if d == days(k)
        pAllStates(k,:) = pNext;
        for e=1:numOfExperiments
            pModel(k,:,e) = pNext * EMatrix(:,:,k,e);
        end
        if derCompFlag
            for i=1:numOfParams
                for e=1:numOfExperiments
                    derpModel(k,:,i,e) = derpNext(:,:,i) * EMatrix(:,:,k,e);
                end
                for j=1:numOfParams
                    for e=1:numOfExperiments
                        secDerpModel(k,:,i,j,e) = secondDerpNext(:,:,i,j) * EMatrix(:,:,k,e);
                    end
                end    
            end  
        end    
        k = k+1;     
    end
    pPrev = pNext;
    %update derivatives
    if derCompFlag
        for i=1:numOfParams
                derpPrev(:,:,i) = derpNext(:,:,i);
                %update second derivatives
                for j=1:numOfParams
                    secondDerpPrev(:,:,i,j) = secondDerpNext(:,:,i,j);
                end
        end
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
        %reshape dimensions of der matrix for param i
        derpMTempI = reshape(derpModel(:,:,i,:), [size(derpModel(:,:,i,:), 1), size(derpModel(:,:,i,:), 2), size(derpModel(:,:,i,:), 4)]); 
        derlogLik(i) = nansum(nansum(nansum(obsMatrix .* derpMTempI ./ pModel)));
        %second derivatives of the log-likelihood
        for j=1:numOfParams
            %reshape dimensions of der and secDer matrices for params i,j
            derpMTempJ = reshape(derpModel(:,:,j,:), [size(derpModel(:,:,j,:), 1), size(derpModel(:,:,j,:), 2), size(derpModel(:,:,j,:), 4)]); 
            secDerpMTempIJ = reshape(secDerpModel(:,:,i,j,:), [size(secDerpModel(:,:,i,j,:), 1), size(secDerpModel(:,:,i,j,:), 2), size(secDerpModel(:,:,i,j,:), 5)]); 
            secDerlogLik(i,j) = nansum(nansum(nansum( obsMatrix .* (secDerpMTempIJ .* pModel - derpMTempI .* derpMTempJ) ./ pModel.^2 )));
        end    
    end
end

%transform to a minimization problem
nL = -logLik;
der = -derlogLik;
%update the hessian in the global variable
Hessian = -secDerlogLik;


end

