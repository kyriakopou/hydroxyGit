function p = oxbisolve(pbi, pox, EBis, EOx)
%OXBISOLVE Solve equation system to get initial distribution from two data
%sets; A * p = [pbi pox 1]

b = [pbi, pox, 1];
A = [EBis, EOx, ones(9,1)];

%use lsqnonneg because the data does not come
%from the same cells and will have inaccuracies
options = optimset('Tolx', 1e-4);
p = lsqnonneg(A', b', options);

end

