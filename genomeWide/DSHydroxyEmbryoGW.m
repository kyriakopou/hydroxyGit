function [nL, der, pBis, pOx, pBisModel, pOxModel, pAllStates] = DSHydroxyEmbryoGW(x, p0, days, obsBis, obsOx, ...
    EBis, EOx, process, derCompFlag)
%EIAP_DSHYDROXY summary of this function goes here
%   Detailed explanation goes here

% x(1): maintenance b0
% x(2): maintenance b1
% x(3): de-Novo b0
% x(4): de-Novo b1
% x(5): hydroxylation b0
% x(6): hydroxylation b1 
% x(7): prob h is recognized as u

global hessian

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

numOfDataPoints = size(days, 1);
numObsStates = size(obsBis, 2);
numOfStates = 9;
numOfParams = size(x, 2);
numOfV = 7;

%modify process vector to get entries 
%for all dataPoints
processAll = cell(days(end), 1);

for i=2:numOfDataPoints
    processAll{days(i)} = process(i);
end
for i=days(end)-1:-1:1
    if (isempty(processAll{i}))
        processAll{i} = processAll{i+1};
    end
end    
% Create constant processes matrices
D = cellDivision(s);


%Bisulfite normalized Observations for all time points 
pBis = obsBis ./ repmat(sum(obsBis, 2), 1, numObsStates);
%Oxidative Bisulfite normalized Observations for all time points
pOx = obsOx ./ repmat(sum(obsOx, 2), 1, numObsStates);


%------Estimate the transient probabilites for each observable day
%Bisulfite - OxBisulfite
%the estimation for p0 are the normalized frequencies of the observations
pBisModel(1,:) = p0 * EBis(:,:,1);
pOxModel(1,:) = p0 * EOx(:,:,1);
pAllStates = zeros(numOfDataPoints, numOfStates);
%the derivatives matrices
derpBisModel = zeros(numOfDataPoints, numObsStates, numOfParams);
derpOxModel = zeros(numOfDataPoints, numObsStates, numOfParams);
%second partial derivatives matrices
secondDerpBisModel = zeros(numOfDataPoints, numObsStates, numOfParams, numOfParams);
secondDerpOxModel = zeros(numOfDataPoints, numObsStates, numOfParams, numOfParams);
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
for d=1:days(end)
    %mu2 and matrices M, derM change at each day bec. mu2 changes
    mu1 = x(1) + x(2) * d;
    mu2 = x(3) + x(4) * d;
    h = x(5) + x(6) * d;
    [M, derM, secondDerM] = methylation(mu1, mu2, p_u, d, derCompFlag);
    [Md, ~, ~] = deNovo(mu2, d, derCompFlag);
    [H, derH, secondDerH] = hydroxylation(h, d, derCompFlag);
    %Compute P Matrix for day d
    if (strcmp(processAll{d}, 'rep'))       
        P = D * M * H;
    elseif (strcmp(processAll{d}, 'no-rep'))
        P = Md * H; 
    end    
    %compute the derivatives of P matrix wrt to every parameter x_i
    if derCompFlag
        for i=1:numOfV
            if (strcmp(processAll{d}, 'rep'))
                derP(:,:,i) = D * (derM(:,:,i)*H + M*derH(:,:,i));
                %compute the second partial derivative (i,j) of P matrix 
                for j=1:numOfV
                    secondDerP(:,:,i,j) = D * (secondDerM(:,:,i,j)*H + derM(:,:,i)*derH(:,:,j) + ... 
                    derM(:,:,j)*derH(:,:,i) + M*secondDerH(:,:,i,j));     
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
   
    %if d = observation day compute the prob. of each of the observations
    %for Bis and OxBis
    if d == days(k)
        pAllStates(k,:) = pNext;
        pBisModel(k,:) = pNext * EBis(:,:,k);
        pOxModel(k,:) = pNext * EOx(:,:,k);
        if derCompFlag
            for i=1:numOfParams
                derpBisModel(k,:,i) = derpNext(:,:,i) * EBis(:,:,k);
                derpOxModel(k,:,i) = derpNext(:,:,i) * EOx(:,:,k);
                for j=1:numOfParams
                    secondDerpBisModel(k,:,i,j) = secondDerpNext(:,:,i,j) * EBis(:,:,k);
                    secondDerpOxModel(k,:,i,j) = secondDerpNext(:,:,i,j) * EOx(:,:,k);
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


%compute likelihood summing over all states over all time points
logLik = nansum(nansum(obsBis .* log(pBisModel) + obsOx .* log(pOxModel)));
%First-second derivatives of the log-likelihood
derlogLik = zeros(1, numOfParams);
secondDerlogLik = zeros(numOfParams, numOfParams);
if derCompFlag
    for i=1:numOfParams
        derlogLik(i) = nansum(nansum(obsBis .* derpBisModel(:,:,i) ./ pBisModel + obsOx .* derpOxModel(:,:,i) ./ pOxModel));
        %second derivatives of the log-likelihood
        for j=1:numOfParams
            secondDerlogLik(i,j) = nansum(nansum( obsBis .* (secondDerpBisModel(:,:,i,j) .* ... 
                pBisModel - derpBisModel(:,:,i) .* derpBisModel(:,:,j)) ./ pBisModel.^2 ...
            + obsOx .* (secondDerpOxModel(:,:,i,j) .* pOxModel - derpOxModel(:,:,i) .* derpOxModel(:,:,j))  ./ pOxModel.^2));
        end    
    end
end

%transform to a minimization problem
nL = -logLik;
der = -derlogLik;
%update the hessian in the global variable
hessian = -secondDerlogLik;


end

