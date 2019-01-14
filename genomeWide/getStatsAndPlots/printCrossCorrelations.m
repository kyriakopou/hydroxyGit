function printCrossCorrelations(M)

%scale the 
Mnums = M(:,2:end);
Mnums_scaled = (Mnums - nanmean(Mnums, 1)) ./ nanstd(Mnums, 1);

Mnums_scaled(isnan(Mnums_scaled)) = 0;
[cross,lags] = xcorr(Mnums_scaled(:,1),'coeff', 10);

numOfCols = size(Mnums_scaled, 2);

for row = 1:numOfCols
    for col = 1:numOfCols
        nm = numOfCols*(row-1)+col;
        subplot(numOfCols,numOfCols,nm)
        cross((size(cross)+1) / 2) = [];
        lags((size(cross)+1) / 2) = [];
        stem(lags,cross(:,nm))
        title(sprintf('c_{%d%d}',row,col))
        ylim([0 1])
    end
end

end