function analyzeOutput(outputMatFile)


if strcmp(outputMatFile(end-5:end-4), 'KO')
 %get the daysLabels from input file
    daysLabels = {'day0', 'day4', 'day7'};
    daysLabelsShort = {'d0', 'd4', 'd7'};
else
    daysLabels = {'day0', 'day3', 'day6'};
    daysLabelsShort = {'d0', 'd3', 'd6'};
end
%filter and reorder the chromosomes
T_out = loadTable(outputMatFile);
% T_out = filterAndReorderChromosomes(T_out);

%getPlots from output file
%%
%avg efficiencies
figure();
plotAvgEfficiencies(T_out, daysLabelsShort);

%xlabel
xticklabels(daysLabels);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);

leg = legend('\mu_m', '\mu_d', '\eta');
leg.FontSize = 15;
ylabel('Efficiency', 'FontSize', 18);

if strcmp(outputMatFile(end-5:end-4), 'KO')
    print('~/Desktop/eff_MLE_KO', '-dpdf');
else
    print('~/Desktop/eff_MLE', '-dpdf');
end    
%%
%plot hidden states over time
figure();
plotHiddenStatesOverTime(T_out, daysLabelsShort);

%xlabel
xticklabels(daysLabels);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 15);

ylabel('Level', 'FontSize', 18);
ylim([0, 1]);
leg = legend({'mm', 'toth', 'hemi', 'uu'}, 'Location', 'north');
leg.FontSize = 15;

if strcmp(outputMatFile(end-5:end-4), 'KO')
    print('~/Desktop/hidden_MLE_KO', '-dpdf');
else
    print('~/Desktop/hidden_MLE_WT', '-dpdf');
end    

%%
%eficiencies per chromosome
uniqueChromosomes = unique(T_out.Chromosome, 'stable');
numOfChromosomes = length(uniqueChromosomes);

fig = figure('Name', 'Efficiencies per Chr',...
    'units', 'centimeters', ...
    'Position', [0, 0, 30, 15]);
for i=1:numOfChromosomes
    subplot(ceil(numOfChromosomes/7), 7, i);   
    rowsChr = strcmp(T_out.Chromosome, uniqueChromosomes{i});
    plotAvgEfficiencies(T_out(rowsChr,:), daysLabelsShort);

    title(uniqueChromosomes{i});
%     ylabel('Efficiency');
    ylim([0,1]);
    
end

if strcmp(outputMatFile(end-5:end-4), 'KO')
    print('-bestfit', '~/Desktop/effPerChr_MLE_KO', '-dpdf');
else
    print('-bestfit', '~/Desktop/effPerChr_MLE_WT', '-dpdf');
end    



%%
%hidden states levels per chromosome 
fig = figure('Name', 'Levels per Chr',...
    'units', 'centimeters', ...
    'Position', [0, 0, 30, 15]);

for i=1:numOfChromosomes
    subplot(ceil(numOfChromosomes/7), 7, i); 
    rowsChr = strcmp(T_out.Chromosome, uniqueChromosomes{i});
    plotHiddenStatesOverTime(T_out(rowsChr,:), daysLabelsShort);

    title(uniqueChromosomes{i});
    ylim([0,1]);
%     ylabel('Level');
%     legend('mm', 'toth', 'um-mu', 'uu');

end

if strcmp(outputMatFile(end-5:end-4), 'KO')
    print('-bestfit', '~/Desktop/hiddenPerChr_MLE_KO', '-dpdf');
else
    print('-bestfit', '~/Desktop/hiddenPerChr_MLE_WT', '-dpdf');
end    

close all;


end