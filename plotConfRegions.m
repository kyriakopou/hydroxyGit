function plotConfRegions(params, Cov)

h = figure();
for i = 1:2:size(params,2)-2
    for j = i+2:2:size(params,2)
        paramsEllipse = [params(i), params(j)];
        CovEllipse = [Cov(i,i), Cov(i,j); Cov(j,i), Cov(j,j)];
        confidenceEllipse(h, paramsEllipse, CovEllipse);
    end   
end

end