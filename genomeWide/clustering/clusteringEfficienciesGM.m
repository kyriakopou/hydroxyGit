function T_clusters = clusteringEfficienciesGM(arg)
%CLUSTERINGEFFICIENCIES Summary of this function goes here
%   Detailed explanation goes here
%The function reads from the WT or TET results files and clusters every CpG
%wrt the enzymes efficiencies, methylation levels, or both

%read from the mat file
filePath = '~/Desktop/epigenetics/wholeGenome/BI_MLE_Results/WT/subtablesMatFiles/';
% fileName = strcat(filePath, 'T_resultsBI_16');
fileName = strcat(filePath, 'T_resultsBI_1_20_XY');


%load output file 
T = load(fileName); 
T_final = T.T_resultsBI_WG;
% T_final = T.T_resultsMLE_WG;
% T_final = T.T_resultsBIChr{16};

%load input file (to read the numOfObs)
inputFile = 'T_input.mat';
T = load(inputFile);
T_input = T.T_final; 

firstDayRows = ~isnan(T_final.maint_d0);
secDayRows = ~isnan(T_final.maint_d3);
thirdDayRows = ~isnan(T_final.maint_d6);

oneDayRows = firstDayRows + secDayRows + thirdDayRows == 1;
twoDaysRows = firstDayRows + secDayRows + thirdDayRows == 2;
threeDaysRows = firstDayRows + secDayRows + thirdDayRows == 3;

%u contains all chromosomes that are not in the results T_final
u = unique(T_input.chr);
u = u(21:end-2);

%delete the corr. rows
f = @(x) strcmp(T_input.chr, x);
del = cellfun(f, u, 'UniformOutput', false);
del = horzcat(del{:});
del = any(del, 2);
T_input(del,:) = [];

%keep only three days CpGs from T_BI
T_final = T_final(threeDaysRows,:);
T_input = T_input(threeDaysRows,:);
%get the number of observations for each CpG
numOfObsPerCpG = sum(T_input{:,3:end}, 2);

clear T;
eff = zeros(size(T_final, 1), 2);

numOfChromosomes = size(T_final, 1);
%read the efficiencies at two timePoints from OutPutFile
%and compute the coefficients (b0, b1) for each of them
eff(:,1) = T_final{:,17};
eff(:,2) = (T_final{:,20} - T_final{:,17}) / 3;

% eff(:,3) = T_final{:,18};
% eff(:,4) = (T_final{:,21} - T_final{:,18}) / 3;
% 
% eff(:,5) = T_final{:,19};
% eff(:,6) = (T_final{:,24} - T_final{:,19}) / 3;
% 
% eff(:,7) = T_final{:,26};

%get levels from T_BI
levels = T_final{:,5:16};

if strcmp(arg, 'both')
    X = [eff, levels];
elseif strcmp(arg, 'eff') 
    X = eff;
else
    X = levels;
end

%Apply Gaussian Mixture Models
maxNumOfClusters = 7;
GMModels = cell(1, maxNumOfClusters);
AIC = zeros(1, maxNumOfClusters);
for k = 1:maxNumOfClusters
    options = statset('Display', 'final', 'MaxIter', 1000);
    GMModels{k} = fitgmdist(X, k, 'Options',options);
    AIC(k)= GMModels{k}.AIC;    
end

[minAIC , numOfClusters] = min(AIC);
numOfClusters

bestModel = GMModels{numOfClusters}
idx = cluster(bestModel, X);
figure();

for k=1:numOfClusters
        rowsCluster = idx == k;
        
        subplot(2, ceil(numOfClusters / 2), k);
        
        if strcmp(arg, 'eff')
            %define colors
            rgbColor1 = [0.8500   0.3250    0.0980;
                        0         0.4470    0.7410;
                        0.9290    0.6940    0.1250];
                    
%             rgbColor1 = [0         0.4470    0.7410];
            
            effAllCluster = {[T_final(rowsCluster,:).maint_d0, T_final(rowsCluster,:).maint_d3, T_final(rowsCluster,:).maint_d6]; 
            [T_final(rowsCluster,:).deNovo_d0, T_final(rowsCluster,:).deNovo_d3, T_final(rowsCluster,:).deNovo_d6];
            [T_final(rowsCluster,:).hydroxy_d0, T_final(rowsCluster,:).hydroxy_d3, T_final(rowsCluster,:).hydroxy_d6]};
            
            %ONLY HYDROXY
%             effAllCluster = {[T_final(rowsCluster,:).hydroxy_d0, T_final(rowsCluster,:).hydroxy_d3, T_final(rowsCluster,:).hydroxy_d6]};          
%             aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd3', 'd6'});
            
            %ONLY DE-NOVO
%             effAllCluster = {[T_final(rowsCluster,:).deNovo_d0, T_final(rowsCluster,:).deNovo_d3, T_final(rowsCluster,:).deNovo_d6]};

            aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd3', 'd6'});
            
        elseif strcmp(arg, 'levels')
            %define colors
            rgbColor2 = [0.8500    0.3250    0.0980;
                        0.929      0.694     0.125;
                        0.466      0.674     0.1880;
                        0          0.447     0.7410];
                        
            levelsAllCluster = {[T_final(rowsCluster,:).mm_d0, T_final(rowsCluster,:).mm_d3, T_final(rowsCluster,:).mm_d6]; 
            [T_final(rowsCluster,:).toth_d0, T_final(rowsCluster,:).toth_d3, T_final(rowsCluster,:).toth_d6];
            [T_final(rowsCluster,:).um_d0, T_final(rowsCluster,:).um_d3, T_final(rowsCluster,:).um_d6]
            [T_final(rowsCluster,:).uu_d0, T_final(rowsCluster,:).uu_d3, T_final(rowsCluster,:).uu_d6]};
        
            aboxplot(levelsAllCluster, 'Colormap', rgbColor2, 'labels', {'d0', 'd3', 'd6'});
        else
            %define colors
            rgbColor1 = [0.8500   0.3250    0.0980;
                        0         0.4470    0.7410;
                        0.9290    0.6940    0.1250];
            
            effAllCluster = {[T_final(rowsCluster,:).maint_d0, T_final(rowsCluster,:).maint_d3, T_final(rowsCluster,:).maint_d6]; 
            [T_final(rowsCluster,:).deNovo_d0, T_final(rowsCluster,:).deNovo_d3, T_final(rowsCluster,:).deNovo_d6];
            [T_final(rowsCluster,:).hydroxy_d0, T_final(rowsCluster,:).hydroxy_d3, T_final(rowsCluster,:).hydroxy_d6]};
        
            aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd3', 'd6'});
            
        end        
              
        ylim([0,1]);
              
        %display numOfCpGs in each cluster
        numOfCpGsInCluster = sum(rowsCluster, 1);
        avgSamplesPerCpG = sum(numOfObsPerCpG(rowsCluster)) / numOfCpGsInCluster;
        fprintf('There are %d CpGs in cluster %d with avg. sample number %.2f \n', numOfCpGsInCluster, k, avgSamplesPerCpG);
        
end

if strcmp(arg, 'both')
figure();
    for k=1:numOfClusters
        
        rowsCluster = idx == k;

        %define colors
        rgbColor2 = [0.8500   0.3250    0.0980;
                    0         0.4470    0.7410;
                    0.466     0.674     0.1880;
                    0         0.447     0.7410];


        levelsAllCluster = {[T_final(rowsCluster,:).mm_d0, T_final(rowsCluster,:).mm_d3, T_final(rowsCluster,:).mm_d6]; 
        [T_final(rowsCluster,:).toth_d0, T_final(rowsCluster,:).toth_d3, T_final(rowsCluster,:).toth_d6];
        [T_final(rowsCluster,:).um_d0, T_final(rowsCluster,:).um_d3, T_final(rowsCluster,:).um_d6]
        [T_final(rowsCluster,:).uu_d0, T_final(rowsCluster,:).uu_d3, T_final(rowsCluster,:).uu_d6]};

        aboxplot(levelsAllCluster, 'Colormap', rgbColor2, 'labels', {'d0', 'd3', 'd6'});
        
        
    end    
end

%empty methylation column
Feature = zeros(size(T_final, 1), 1);
%store cluster for every CpG in T_out and copy this to .mat, .txt, .igv file 
varNames = {'Chromosome', 'Start', 'End', 'Methylation', 'Cluster', 'numOfDataPoints'};
T_clusters = table(T_final.Chromosome, T_final.Start, T_final.End, Feature, idx, numOfObsPerCpG, 'VariableNames', varNames);
save(strcat(filePath, 'T_clusteringBI'), 'T_clusters');

%write results of MLE table to a .txt file
txt_FilePath = strcat(filePath, 'resultsClusteringBI.txt');
writetable(T_clusters, txt_FilePath, 'Delimiter', '\t')
% %change the MLE file's extension from txt to .igv
% igv_FilePath = strcat(filePath, 'resultsClusteringBI.igv');
% movefile(txt_FilePath, igv_FilePath);

end

