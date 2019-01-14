function analyzeInput(inputMatFile)


if strcmp(inputMatFile(end-5:end-4), 'KO')
 %get the daysLabels from input file
    daysLabels = {'day0', 'day4', 'day7'};
    daysLabelsPerChr = {'d0', 'd4', 'd7'};
else
    daysLabels = {'day0', 'day3', 'day6'};
    daysLabelsPerChr = {'d0', 'd3', 'd6'};
end

%get the table of the inputfile
T_in = loadTable(inputMatFile);

%get statistics of the data in T_final
[oneDayRows, twoDaysRows, threeDaysRows] = getRowsOfTable(T_in);

numOfCpGsWithOneTimePoints = nnz(oneDayRows);
numOfCpGsWithTwoTimePoints = nnz(twoDaysRows);
numOfCpGsWithThreeTimePoints = nnz(threeDaysRows);
numOfCpGsPerNumOfTimePoints = [numOfCpGsWithOneTimePoints, numOfCpGsWithTwoTimePoints, numOfCpGsWithThreeTimePoints];

%plot the numOfCpGs with 1, 2, 3 timepoints
fig = figure();
plotNumOfCpGsDetailed([numOfCpGsPerNumOfTimePoints]);
title('Number of CpGs measured', 'FontSize', 20);

%xaxis fontsize
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
%yaxis
ylabel('Number of CpGs', 'FontSize', 18);
%print pdf
print('~/Desktop/numObsDays', '-dpdf');

%%
%extend T_in to contain the total number of samples per day (BS and oxBS) for each CpG 
T_in = createTotalObsPerDayTable(T_in);
%print it out
% printSamplesPerChrTable(T_in);
%plot avg number of samples per day gw
avgObs = avgNumOfSamplesPerDay(T_in);

figure();
axes('XTick', 1:3, 'XTickLabel', daysLabels);
% plotAvgNumOfSamplesPerDay(avgObs);
plotNumOfSamplesPerDay(T_in, daysLabels);

title('Number of samples per day', 'FontSize', 20);

%xaxis
xticklabels(daysLabels);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);
%yaxis
ylabel('Number of samples', 'FontSize', 18);
leg = legend('BS', 'oxBS', 'Location', 'northeast');
leg.FontSize = 15;

% set(findall(gcf,'-property','FontSize'),'FontSize', 15)
print('~/Desktop/numOfSamples', '-dpdf');


%%
%filter and reorder the chromosomes
T_in = filterAndReorderChromosomes(T_in);

uniqueChromosomes = unique(T_in.Chromosome, 'stable');
numOfChromosomes = size(uniqueChromosomes, 1);

%plot numOfCpGs with 1, 2, 3 timepoints for all chromosomes
fig = figure('Name', 'Number Of CpGs measured',...
    'units', 'centimeters', 'Position', [0, 0, 30, 15]);

for i=1:numOfChromosomes
    
    rowsChr = strcmp(T_in.Chromosome, uniqueChromosomes{i});
    [oneDayRows, twoDaysRows, threeDaysRows] = getRowsOfTable(T_in(rowsChr,:));

    numOfCpGsWithOneTimePoints = nnz(oneDayRows);
    numOfCpGsWithTwoTimePoints = nnz(twoDaysRows);
    numOfCpGsWithThreeTimePoints = nnz(threeDaysRows);
    numOfCpGsPerNumOfTimePoints = [numOfCpGsWithOneTimePoints, numOfCpGsWithTwoTimePoints, numOfCpGsWithThreeTimePoints];

    subplot(ceil(numOfChromosomes/7), 7, i);
    plotNumOfCpGsDetailed(numOfCpGsPerNumOfTimePoints);
    title({uniqueChromosomes{i}; ''});
   
end
print('-bestfit', '~/Desktop/numObsDaysPerChr', '-dpdf');

%%
%plot avg number of samples per day and chromosome
fig = figure('Name', 'NumOfSamples per day for all Chrs',...
    'units', 'centimeters', 'Position', [0, 0, 30, 15]);

for i=1:numOfChromosomes
    
    subplot(ceil(numOfChromosomes/7), 7, i);
    rowsChr = strcmp(T_in.Chromosome, uniqueChromosomes{i});
%     plotAvgNumOfSamplesPerDay(avgNumOfSamplesPerDay(T_in(rowsChr,:)));
    plotNumOfSamplesPerDay(T_in(rowsChr,:), daysLabelsPerChr); 
    xticklabels(daysLabelsPerChr);
    title(uniqueChromosomes{i});
    
end
print('-bestfit', '~/Desktop/numOfSamplesPerChr', '-dpdf');

close all;

end