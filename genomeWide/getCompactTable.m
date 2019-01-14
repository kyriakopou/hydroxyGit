function T_compact = getCompactTable(matDirName)

extTableMatFileName = strcat(matDirName, '/T_BI_WG.mat');
T = load(extTableMatFileName);
%change this to T_all_mle if desired
T_all = T.T_all_bi; 
clear T;

%extract observation days from the variable names
days = str2num([T_all.Properties.VariableNames{5}(end); T_all.Properties.VariableNames{14}(end); T_all.Properties.VariableNames{23}(end)])';
numOfDays = length(days);
numOfColsInCompactMatrix = numOfDays*4 + (numOfDays*3+1)*2;
compArray = nan(size(T_all, 1), numOfColsInCompactMatrix);

%copy the parts of the table u need to huge matrices
pAllStatesColumns = T_all{:, 5:31};
paramsColumns = T_all{:, 32:38};
CovUpTriangleColumns = T_all{:,39:66};

% for i=1:size(T_all, 1)
% profile on;
tic;
parfor i=2946147:size(T_all, 1)
    %get the distribution of the states
    pAllStates = pAllStatesColumns(i,:);
    pAllStates = reshape(pAllStates, [9, numOfDays])';
    
    %get summarized distribution
    uuCol = pAllStates(:,1);
    umCol = sum(pAllStates(:,2:3), 2);
    tothCol = sum(pAllStates(:,4:5), 2) + sum(pAllStates(:,6:7), 2) + pAllStates(:,9);
    mmCol = pAllStates(:,8);
    
    sumDist = [mmCol, tothCol, umCol, uuCol];
    
    %extract params and Cov matrix
    x = paramsColumns(i,:);
    
    CovUpTriangle = CovUpTriangleColumns(i,:);
    numOfParams = numel(x);
    Cov = zeros(numOfParams, numOfParams);
    Cov(triu(true(numOfParams))) = CovUpTriangle;
    Cov = Cov + tril(Cov', -1);   
    %add the missing entries to the cov needed for getEfficiencies function
    CovEff = zeros(11,11);
    CovEff(1:6, 1:6) = Cov(1:6, 1:6);

    %call getEfficiencies to get the values of the rates and their covs at days 
    params = [x(1:6) zeros(1,4) x(7)];
    [Eff, varEff] = getEfficiencies(params, 1, days, [], CovEff);
%     [Eff, varEff] = getEfficiencies(params, 1, days, [], CovPlot);
    Eff = Eff(:,1:3);
    varEff = varEff(:,1:3);
    
    %create array with compact results
    compArray(i,:) = [reshape(sumDist', [1 numOfDays * 4]), reshape(Eff', [1, numOfDays * 3]), x(end), reshape(sqrt(varEff)', [1, numOfDays * 3]), Cov(7,7)];
    
end

toc;
headersNumer = {'mm_d0', 'toth_d0', 'um_d0' 'uu_d0', 'mm_d3', 'toth_d3', 'um_d3', 'uu_d3', 'mm_d6', 'toth_d6', 'um_d6', ...
    'uu_d6', 'maint_d0', 'deNovo_d0', 'hydroxy_d0', 'maint_d3', 'deNovo_d3', 'hydroxy_d3', 'maint_d6', 'deNovo_d6', ...
    'hydroxy_d6', 'pRecogn', 'maint_d0_sd', 'deNovo_d0_sd', 'hydroxy_d0_sd', 'maint_d3_sd' , 'deNovo_d3_sd', 'hydroxy_d3_sd', ...
    'maint_d6_sd', 'deNovo_d6_sd', 'hydroxy_d6_sd', 'pRecogn_sd'};

% profile viewer;
% profile off;

T_compact = [T_all(:,1:4), array2table(compArray, 'VariableNames', headersNumer)];
save(strcat(matDirName, '/T_compact_BI_WG_.mat'), 'T_compact');

end