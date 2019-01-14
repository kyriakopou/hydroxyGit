function [CpGs, idx] = clusteringEfficienciesTetKO(arg)
%CLUSTERINGEFFICIENCIES Summary of this function goes here
%   Detailed explanation goes here
%The function reads from the WT or TET results files and clusters every CpG
%wrt the enzymes efficiencies, methylation levels, or both

%read from the mat file
filePath = '~/Desktop/epigenetics/wholeGenome/BI_MLE_Results/TET_KO/subtablesMatFiles/';
% fileName = strcat(filePath, 'T_resultsBI_16');
% fileName = strcat(filePath, 'T_resultsBI_1_20_XY');
%TET_KO example
fileName = strcat(filePath, 'T_resultsBI_WG');


%load output file 
T = load(fileName); 
T_final = T.T_resultsBI_WG;
% T_final = T.T_resultsMLE_WG;
% T_final = T.T_resultsBIChr{16};

%load input file (to read the numOfObs)
inputFile = 'T_inputK0.mat';
T = load(inputFile);
T_input = T.T_final; 

firstDayRows = ~isnan(T_final.maint_d0);
secDayRows = ~isnan(T_final.maint_d4);
thirdDayRows = ~isnan(T_final.maint_d7);

oneDayRows = firstDayRows + secDayRows + thirdDayRows == 1;
twoDaysRows = firstDayRows + secDayRows + thirdDayRows == 2;
threeDaysRows = firstDayRows + secDayRows + thirdDayRows == 3;



%keep only three days CpGs from T_BI
T_final = T_final(threeDaysRows,:);
T_input = T_input(threeDaysRows,:);
%get the number of observations for each CpG
numOfObsPerCpG = sum(T_input{:,3:end}, 2);

clear T;
eff = zeros(size(T_final, 1), 2);

numOfChromosomes = size(T_final, 1);
% for i=1:numOfChromosomes
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

%choose k and
%make k-means clustering based on the efficiencies
%     numOfClusters = 6;
%     idx = kmeans(eff, numOfClusters);
%evaluate the clusterings for a specific range of clusters
myFunc = @(X,K)(kmeans(X, K, 'emptyaction', 'singleton',...
    'replicate', 10));
eva = evalclusters(X, myFunc, 'CalinskiHarabasz', 'KList', 2:10);
plot(eva);
numOfClusters = eva.OptimalK;
idx = eva.OptimalY;
% end


figure();
for k=1:numOfClusters
        rowsCluster = idx == k;
        
        subplot(2, ceil(numOfClusters / 2), k);
        
        if strcmp(arg, 'eff')
            %define colors
            rgbColor1 = [0.8500   0.3250    0.0980;
                        0         0.4470    0.7410];
            
            effAllCluster = {[T_final(rowsCluster,:).maint_d0, T_final(rowsCluster,:).maint_d4, T_final(rowsCluster,:).maint_d7]; 
            [T_final(rowsCluster,:).deNovo_d0, T_final(rowsCluster,:).deNovo_d4, T_final(rowsCluster,:).deNovo_d7]};
            
            %ONLY HYDROXY
%             effAllCluster = {[T_final(rowsCluster,:).hydroxy_d0, T_final(rowsCluster,:).hydroxy_d3, T_final(rowsCluster,:).hydroxy_d6]};          
%             aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd3', 'd6'});
            
            %ONLY DE-NOVO
%             effAllCluster = {[T_final(rowsCluster,:).deNovo_d0, T_final(rowsCluster,:).deNovo_d3, T_final(rowsCluster,:).deNovo_d6]};

            aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd4', 'd7'});
            
        elseif strcmp(arg, 'levels')
            %define colors
            rgbColor2 = [0.8500    0.3250    0.0980;
                        0.466      0.674     0.1880;
                        0          0.447     0.7410];
                        
            levelsAllCluster = {[T_final(rowsCluster,:).mm_d0, T_final(rowsCluster,:).mm_d4, T_final(rowsCluster,:).mm_d7]; 
            [T_final(rowsCluster,:).um_d0, T_final(rowsCluster,:).um_d4, T_final(rowsCluster,:).um_d7]
            [T_final(rowsCluster,:).uu_d0, T_final(rowsCluster,:).uu_d4, T_final(rowsCluster,:).uu_d7]};
        
            aboxplot(levelsAllCluster, 'Colormap', rgbColor2, 'labels', {'d0', 'd4', 'd7'});
        else
            %define colors
            rgbColor1 = [0.8500   0.3250    0.0980;
                        0         0.4470    0.7410];
            
            effAllCluster = {[T_final(rowsCluster,:).maint_d0, T_final(rowsCluster,:).maint_d4, T_final(rowsCluster,:).maint_d7]; 
            [T_final(rowsCluster,:).deNovo_d0, T_final(rowsCluster,:).deNovo_d4, T_final(rowsCluster,:).deNovo_d7]};
        
            aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd4', 'd7'});
            
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
                    0.466     0.674     0.1880;
                    0         0.447     0.7410];


        levelsAllCluster = {[T_final(rowsCluster,:).mm_d0, T_final(rowsCluster,:).mm_d4, T_final(rowsCluster,:).mm_d7]; 
        [T_final(rowsCluster,:).um_d0, T_final(rowsCluster,:).um_d4, T_final(rowsCluster,:).um_d7]
        [T_final(rowsCluster,:).uu_d0, T_final(rowsCluster,:).uu_d4, T_final(rowsCluster,:).uu_d7]};

        aboxplot(levelsAllCluster, 'Colormap', rgbColor2, 'labels', {'d0', 'd4', 'd7'});
        
        
    end    
end

%empty methylation column
Feature = zeros(size(T_final, 1), 1);
%store cluster for every CpG in T_out and copy this to .mat, .txt, .igv file 
varNames = {'Chromosome', 'Start', 'End', 'Methylation', 'Cluster', 'numOfDataPoints'};
T_out = table(T_final.Chromosome, T_final.Start, T_final.End, Feature, idx, numOfObsPerCpG, 'VariableNames', varNames);
save(strcat(filePath, 'T_clusteringBI'), 'T_out');



%write results of MLE table to a .txt file
txt_FilePath = strcat(filePath, 'resultsClusteringBI.txt');
writetable(T_out, txt_FilePath, 'Delimiter', '\t')
% %change the MLE file's extension to .igv
% igv_FilePath = strcat(filePath, 'resultsClusteringBI.igv');
% movefile(txt_FilePath, igv_FilePath);

end
