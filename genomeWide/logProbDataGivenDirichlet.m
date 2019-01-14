function logProb = logProbDataGivenDirichlet(X_bis, X_ox, A_bis, A_ox, numOfDays)

%get the sums of each row of params array
a0_bis = sum(A_bis, 2);
a0_ox = sum(A_ox, 2);
%the same for the samples array
nBis = sum(X_bis, 2);
nOx = sum(X_ox, 2);

%split the sum of logs in two terms 
aTermBis = zeros(numOfDays, 1);
bTermBis = zeros(numOfDays, 1);
aTermOx = zeros(numOfDays, 1);
bTermOx = zeros(numOfDays, 1);

%mean of dirichlet
meanBis = repmat(nBis, 1, 4) .* X_bis ./ repmat(a0_bis, 1, 4);
meanOx = repmat(nOx, 1, 4) .* X_ox ./ repmat(a0_ox, 1, 4);
%approximate mode of dirichlet
apprModeBis = round(meanBis);
apprModeOx = round(meanOx);
%prob of the mode of DirMultin
logPrBisMeanOfDay = zeros(numOfDays, 1);
logPrOxMeanOfDay = zeros(numOfDays, 1);


for k=1:numOfDays

        aTermBis(k,:) = sum(log(1:a0_bis(k)-1));
        bTermBis(k,:) = sum(log(nBis(k)+1:nBis(k)+a0_bis(k)-1));
        logPrBisMeanOfDay(k) = dirichletMultinomialPmf(apprModeBis(k, :), A_bis(k,:));
        aTermOx(k,:) = sum(log(1:a0_ox(k)-1));
        bTermOx(k,:) = sum(log(nOx(k)+1:nOx(k)+a0_ox(k)-1));
        logPrOxMeanOfDay(k) = dirichletMultinomialPmf(apprModeOx(k, :), A_ox(k,:));        
end

temp1MatrixBis = zeros(size(X_bis));
temp2MatrixBis = zeros(size(X_bis));
temp1MatrixOx = zeros(size(X_ox));
temp2MatrixOx = zeros(size(X_ox));

%create the entries of tempMatrices (last two terms of the sum)
for row = 1:size(X_bis,1)
    for col = 1:size(X_bis,2)
        elmBis = X_bis(row, col);
        elmOx = X_ox(row, col);

        temp1MatrixBis(row,col) = sum(log(elmBis+1 : elmBis+A_bis(row,col)-1 )); 
        temp2MatrixBis(row,col) = sum(log(1:A_bis(row, col)-1));
        temp1MatrixOx(row,col) = sum(log(elmOx+1 : elmOx+A_ox(row,col)-1 )); 
        temp2MatrixOx(row,col) = sum(log(1:A_ox(row, col)-1));
    end
end
    
%sum up all the derived terms to get the logProb of the current CpG
%data based on the current model (dirichlet params)
%also weight the final prob of each day with the prob of the mean of DirMult
%if you dont want to weight remove logPrBisMeanOfDay and logPrOxMeanOfDay
probBisDayCpGFromGamma = aTermBis - bTermBis + sum(temp1MatrixBis, 2) - sum(temp2MatrixBis, 2) - logPrBisMeanOfDay; 
probOxDayCpGFromGamma = aTermOx - bTermOx + sum(temp1MatrixOx, 2) - sum(temp2MatrixOx, 2) - logPrOxMeanOfDay;

logProb = sum(probBisDayCpGFromGamma + probOxDayCpGFromGamma, 1);

if isnan(logProb)
    disp(logProb);
end    

end    