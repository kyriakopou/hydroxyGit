function runDataGW(dataFilesPath, chromosomesToRun)
%given a fixed clustering of the CpGs run the parameter estimation procedure
%(either MLE or bayesian inference) for each cluster

tic;

switch nargin 
    case 1
        chromosomesToRun = 1;
end        
        
%read the errorsFile
errorFile = strcat(dataFilesPath, '/errors.txt');
[~, maxObsDaysVector, maxDataLabels, errors_BS, errors_oxBS, process] = readErrors(errorFile);

%get the T_final table
T_final = preprocessData(dataFilesPath, maxObsDaysVector);

% T_final.chr = strcat('chr', T_final.chr);
%split the data in different chromosomes
[uniqueChromosomes, indFirstCpG] = findChromosomes(T_final);
numOfChromosomes = size(uniqueChromosomes, 1);

% numOfChromosomes = 1;  %without the strings Chromosomes (X,Y included)
windowSize = 1;
numOfBinsPerChr = zeros(numOfChromosomes, 1);
totNumOfCpGs = size(T_final, 1);

headersStrings = {'Chromosome', 'Start', 'End', 'Methylation'}; 
 
headersNumer = {'uu_d0', 'um_d0', 'mu_d0' 'uh_d0', 'hu_d0', 'hm_d0', 'mh_d0', 'mm_d0', 'hh_d0', 'uu_d3', 'um_d3', ...
    'mu_d3' 'uh_d3', 'hu_d3', 'hm_d3', 'mh_d3', 'mm_d3', 'hh_d3', 'uu_d6', 'um_d6', 'mu_d6' 'uh_d6', 'hu_d6', 'hm_d6', ...
    'mh_d6', 'mm_d6', 'hh_d6', 'maint_b0', 'maint_b1', 'deNovo_b0', 'deNovo_b1', 'hydroxy_b0', 'hydroxy_b1', 'pRecogn', ...
    'cov_11', 'cov_12', 'cov_22', 'cov_13', 'cov_23', 'cov_33', 'cov_14', 'cov_24', 'cov_34', 'cov_44', 'cov_15', ... 
    'cov_25', 'cov_35', 'cov_45', 'cov_55', 'cov_16',  'cov_26', 'cov_36', 'cov_46', 'cov_56', 'cov_66', 'cov_17', ...
    'cov_27','cov_37', 'cov_47', 'cov_57', 'cov_67', 'cov_77', 'KLBis', 'KLOx', 'hyperVol', 'accRate'};


 
%%--FOR THE OUTPUT--
%(1st path is for WT experiments, 2nd path is for tetKO experiment)
% USE /local/charis/genWideResults/winSize_ and /local/charis/genWideResultsTetKO/winSize_ for alma machine
resultsPath = strcat('~/Desktop/genWideResults/winSize_', num2str(windowSize));
% resultsPath = strcat('~/Desktop/genWideResultsTetKO/winSize_', num2str(windowSize));

%create needed directories if they dont exist
mlePath = strcat(resultsPath, '/allChromosomes/MLE');
biPath = strcat(resultsPath, '/allChromosomes/BI');
if ~(exist(mlePath, 'dir') == 7)
    mkdir(mlePath);
end
if ~(exist(biPath, 'dir') == 7)
    mkdir(biPath);
end

%initialize chromosome tables
T_resultsMLEChr = cell(numOfChromosomes, 1);
T_resultsBIChr = cell(numOfChromosomes, 1); 

% numOfCpGsPerChr = getNumOfCpGsPerChromosome(indFirstCpG, totNumOfCpGs);

% for k=1:numOfChromosomes
for k=chromosomesToRun
    
%     profile on;
  
    %get the subtable that refers to uniqueChromosome(k)
    if k < numOfChromosomes
        T_chr = T_final(indFirstCpG(k):indFirstCpG(k+1)-1, :);
    %last chromosome     
    else
        T_chr = T_final(indFirstCpG(k):totNumOfCpGs, :);
    end
    
    %ONLY WHEN RUNNING ONE CHROMOSOME TO SAVE MEMORY
    %FOR THE REST DELETE ONLY THE USED CHROMOSOMES
%     clear T_final;
    chr = T_chr{1,1}{1};
    fprintf('================================chr%s START ====================================\n\n', chr);
    
    %aggregate the data wrt to some windowSize or for given regions-splits
    %for a given clustering of the data C(:,1) is the bin, C(:,2) the CpGs
    %inside this bin --THIS FUNCTION IS A BIT SLOW-MAYBE OPTIMIZE IT
    [obsBis, obsOx, C] = aggregateData(T_chr, windowSize);
    
    %--HAS TO BE CHANGED FOR MORE THAN ONE CpGs in a bin--
    %CpG region(bin) and CpG that contains
    %copy the first and the last CpG of the estimated region for the output to IGV file
    Start = cell2mat(C(:,2))-1;
    End = cell2mat(C(:,2))+1;
    
    numOfBinsPerChr(k) = size(C, 1);
      
    %rearrange obsBis, obsOx columns (if he changes the order of the col remove this)
    obsBis = obsBis(:, [4 3 2 1], :);
    obsOx = obsOx(:, [4 3 2 1], :);
    
    %create a directory for that chromosome if it does not exist
    chrDir = strcat(resultsPath, '/chr', uniqueChromosomes{k});
    if ~(exist(chrDir, 'dir') == 7)
        mkdir(chrDir);
    end
    
    MLEChrFileName = strcat(mlePath, '/T_resultsMLE_Chr_', mat2str(chromosomesToRun), '.mat');
    BIChrFileName = strcat(biPath, '/T_resultsBI_Chr_', mat2str(chromosomesToRun), '.mat');
    
    %run the estimation for all windows (splits)
    %and store the results in output table T_results
    %here we initialize all the variables (columns) of T_results  
    %chromosome array numer results. we store the intermediate results
    %of each split of each chromosome in numResultsMLE-BIArray in case
    %something goes wrong
    Feature = zeros(numOfBinsPerChr(k), 1);
    
    %split in sets of numOfCpGsPerSplit CpGs (bins) if more than 100 bins
    maxNumOfBinsPerSplit = 500;
    numOfSplits = ceil(numOfBinsPerChr(k) / maxNumOfBinsPerSplit);
    
    %read numResBIArray.mat file of the chr and continue the computations from there
    %in case it is not empty
    initialSplit = 1;
    if exist(BIChrFileName, 'file')
        %load old results and copy to numResultsArrays
        S = load(BIChrFileName);
        numResultsBIArray = S.T_resultsBIChr{chromosomesToRun}{:,5:end};
        S = load(MLEChrFileName);
        numResultsMLEArray = S.T_resultsMLEChr{chromosomesToRun}{:,5:end};
        
        %find the first empty row (CpG) to begin from there the computations
        firstEmptyRow = find(all(isnan(numResultsBIArray), 2), 1);

        if ~isempty(firstEmptyRow)
            initialSplit = ceil(firstEmptyRow / maxNumOfBinsPerSplit);
        else
            initialSplit = numOfSplits+1;
        end
    else
        %initialize numResultsArrays to nan
        numResultsMLEArray = nan(numOfBinsPerChr(k), size(headersNumer, 2));
        numResultsBIArray = nan(numOfBinsPerChr(k), size(headersNumer, 2));
    end
    
      
    for s = initialSplit:numOfSplits
              
        %FOR PARALLELIZED VERSION UNCOMMENT THE FOLLOWING LINES
        if isempty(gcp('nocreate'))
            %change the job storage location for each chromosome (each machine
            %in Deep Cluster) to avoid conflicts
            pc = parcluster;
            uniqueJobPath = strcat(pc.JobStorageLocation, '/chr', uniqueChromosomes{k});
            if (~exist(uniqueJobPath, 'dir')) 
                mkdir(uniqueJobPath);
            end
            pc.JobStorageLocation = uniqueJobPath;            
            %parpool(pc, 24);
            parpool(pc);
        end
         
        splitsFirstBin = (s-1)*maxNumOfBinsPerSplit+1;
        splitsLastBin = min(s*maxNumOfBinsPerSplit, numOfBinsPerChr(k));
        numOfBins = splitsLastBin - splitsFirstBin + 1;
        
        %copy only the relevant part (split) of chromosome
        obsBisSplit = obsBis(:,:,splitsFirstBin:splitsLastBin);
        obsOxSplit = obsOx(:,:,splitsFirstBin:splitsLastBin);
        binAndCpGsSplit = C(splitsFirstBin:splitsLastBin,:);
        cpgArray = cell2mat(binAndCpGsSplit(:, 2:end));
        numResultsMLEArraySplit = nan(numOfBins, 66);
        numResultsBIArraySplit = nan(numOfBins, 66);
        
        %run in parallel each CpG of the split          
        parfor i = 1:numOfBins
            
            tStartBin = tic;
                                                
            %call estimateDSHydroxy (parameter estimation using MLE) --
            %give chrDir as the last argument if u want to write the results 
            %to a .txt file
            [params, paramNames, Cov, pAllStates, KLBis, KLOx, hyperVol] = estimateDSHydroxyGW(chr, binAndCpGsSplit(i,:), ...
                maxObsDaysVector, maxDataLabels, obsBisSplit(:,:,i), obsOxSplit(:,:,i), ...
                errors_BS, errors_oxBS, process, chrDir);
            
           
            screenOutput = 0;
            %print CpG info in the screen
            if screenOutput
                %get the string of CpGs inside the bin
                stringCpGs = sprintf('%d_', cpgArray(i,:));
                tEndMLE = toc(tStartBin);
                fprintf('CpG_%s:\t', stringCpGs);
                if all(isnan(params))
                    fprintf('MLE-lev:%4.2f\t', tEndMLE); 
                else
                    fprintf('MLE-eff:%4.2f\t', tEndMLE);
                end
            end
            
            %write detailed results to a matrix
            numResultsMLEArraySplit(i,:) = getDetailedNumerResults(params, Cov, pAllStates, KLBis, KLOx, hyperVol);
            
            %the BI estimation will be executed only in case of 3 obsDays
            if size(find(~all(obsBisSplit(:,:,i) == 0, 2)), 1) == 3
                tStartBI = tic;
                %Run the BI computation %give chrDir as the last argument if u want to write the results 
                %to a .txt file
                [params, Cov, pAllStates, KLBis, KLOx, hyperVol, accRate] = bayesianDSHydroxyGW(binAndCpGsSplit(i,:), ...
                    maxObsDaysVector, maxDataLabels, obsBisSplit(:,:,i), obsOxSplit(:,:,i), errors_BS, errors_oxBS, ...
                    process, pAllStates(1,:), params, Cov, paramNames, chrDir);
                              
                tEndBI = toc(tStartBI);
                if screenOutput
                    fprintf('BI:%4.2f(%.2f)\t', tEndBI, accRate);
                end    
                %write detailed results to a matrix
                numResultsBIArraySplit(i,:) = getDetailedNumerResults(params, Cov, pAllStates, KLBis, KLOx, ...
                    hyperVol, accRate);

            else
                %just copy the MLE estimated levels to BI table
                numResultsBIArraySplit(i,:) = numResultsMLEArraySplit(i,:);
                
                if screenOutput
                    %leave empty horizontal space
                    fprintf('%13s\t', []);   
                end
            end
                      
            if screenOutput
                tFinalBin = toc(tStartBin);
                fprintf('total:%4.2f \n', tFinalBin);
            end
%             parTimes(i) = Par.toc;       
%             fprintf('total: %4.2f sec.\n', parTimes(i).ItStop - parTimes(i).ItStart);
        end
%         stop(parTimes);

        delete(gcp('nocreate'));
        
        %pass to the chromosome array results
        numResultsMLEArray(splitsFirstBin:splitsLastBin, :) = numResultsMLEArraySplit;
        numResultsBIArray(splitsFirstBin:splitsLastBin, :) = numResultsBIArraySplit;
        
         
        lastChrMLEMatFileName = strcat(chrDir, '/numResMLEArray.mat');
        lastChrBIMatFileName = strcat(chrDir, '/numResBIArray.mat');
        %save the results of chromosome so far in .mat format
        save(lastChrMLEMatFileName, 'numResultsMLEArray');
        save(lastChrBIMatFileName, 'numResultsBIArray');
        
        numOfFinishedBins = s*numOfBins;
        if s~= numOfSplits
            fprintf('\n--%d / %d chromosomes and %d / %d CpGs of chr%s have finished--\n', k-1, ...
                numOfChromosomes, numOfFinishedBins, numOfBinsPerChr(k), chr);    
        else
            fprintf('\n--%d / %d chromosomes have finished--\n', k, numOfChromosomes);  
        end
        
    end
    
    %table containing chr, position info
    T_position = table(T_chr{1:numOfBinsPerChr(k),1}, Start(1:numOfBinsPerChr(k)), End(1:numOfBinsPerChr(k)), ... 
        Feature, 'VariableNames', headersStrings);
    %table containing MLE results of this split
    T_numResMLE = array2table(numResultsMLEArray, 'VariableNames', headersNumer); 
    %store in MLE table for the whole chromosome
    T_resultsMLEChr{k} = horzcat(T_position, T_numResMLE);
    %store in BI table for the whole chromosome
    T_numResBI = array2table(numResultsBIArray, 'VariableNames', headersNumer);
    T_resultsBIChr{k} = horzcat(T_position, T_numResBI);
    %save the tables in workspace
    save(MLEChrFileName, 'T_resultsMLEChr');
    save(BIChrFileName, 'T_resultsBIChr');
    

end


%write general information about the experiment
% fileID = fopen('/local/charis/myCode/MATLAB/HydroxyMethylation/genomeWide/results/generalInfo.txt', 'w');
chrGenInfoPath = strcat(chrDir, '/generalInfo_chr', uniqueChromosomes{k}, '.txt');

fileID = fopen(chrGenInfoPath, 'w');
fprintf(fileID, ...
    'The algorithm has been run for chromosomes %s within the genome with \n total number of CpGs equal to %d. ', ...
    uniqueChromosomes{chromosomesToRun}, numOfBinsPerChr(chromosomesToRun));
fprintf(fileID, 'The total running time was %.1f sec.\n\n', toc); 
for k=chromosomesToRun
   fprintf(fileID, '%s:  %d CpG regions \n', uniqueChromosomes{k}, numOfBinsPerChr(k));   
end
fclose(fileID);

%get profiler results
% profile viewer;
% htmlResPath = strcat(resultsPath, '/myProfile'); 
% profsave(profile('info'), htmlResPath);
% profile off;

end
