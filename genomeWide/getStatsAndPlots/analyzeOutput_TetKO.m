function analyzeOutput(outputMatFile)

daysLabels = {'day0', 'day4', 'day7'};

T_out = loadOutputTable(outputMatFile);

%getPlots from output file

%%
%avg efficiencies
figure();
% plotAvgEfficienciesGW(T_out, daysLabels);
plotAvgEfficiencies_TetKO(T_out, daysLabels);

%xlabel
xticklabels(daysLabels);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 12);

leg = legend('\mu_m', '\mu_d', '\eta');
leg.FontSize = 12;
ylabel('Efficiency', 'FontSize', 15);


print('~/Desktop/eff_tet', '-dpdf');
%%
%plot hidden states over time
figure();
% plotHiddenStatesOverTimeGW(T_out, daysLabels);
plotHiddenStatesOverTime_TetKO(T_out, daysLabels);

%xlabel
xticklabels(daysLabels);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 12);

ylabel('Level', 'FontSize', 15);
ylim([0, 1]);
leg = legend({'mm', 'toth', 'hemi', 'uu'}, 'Location', 'north');
leg.FontSize = 12;

print('~/Desktop/hidden_tet', '-dpdf');

%%
%eficiencies per chromosome
[uniqueChromosomes, ~] = findChromosomes(T_out);
numOfChromosomes = size(uniqueChromosomes, 1);
daysLabels = {'d0', 'd4', 'd7'};

fig = figure('Name', 'Efficiencies per Chr',...
    'units', 'centimeters', ...
    'Position', [0, 0, 30, 30]);
for i=1:numOfChromosomes/2
    subplot(ceil(numOfChromosomes/12), 6, i);   
    rowsChr = strcmp(T_out.Chromosome, uniqueChromosomes{i});
%     plotAvgEfficienciesGW(T_out(rowsChr,:), daysLabels);
    plotAvgEfficiencies_TetKO(T_out(rowsChr,:), daysLabels);   

    title(uniqueChromosomes{i});
%     ylabel('Efficiency');
    ylim([0,1]);
    
end
print('-bestfit', '~/Desktop/effPerChr_tet_a', '-dpdf');


fig = figure('Name', 'Efficiencies per Chr',...
    'units', 'centimeters', ...
    'Position', [0, 0, 30, 30]);
for i=1:numOfChromosomes/2
    subplot(ceil(numOfChromosomes/12), 6, i);   
    rowsChr = strcmp(T_out.Chromosome, uniqueChromosomes{i+numOfChromosomes/2});
%     plotAvgEfficienciesGW(T_out(rowsChr,:), daysLabels);   
    plotAvgEfficiencies_TetKO(T_out(rowsChr,:), daysLabels);   

    title(uniqueChromosomes{i+numOfChromosomes/2});
%     ylabel('Efficiency');
    ylim([0,1]);
    
    % legend('mu_m', '\mu_d', '\eta');
end
print('-bestfit', '~/Desktop/effPerChr_tet_b', '-dpdf');


%%
%hidden states levels per chromosome 
fig = figure('Name', 'Levels per Chr',...
    'units', 'centimeters', ...
    'Position', [0, 0, 30, 30]);

for i=1:numOfChromosomes/2
    subplot(ceil(numOfChromosomes/12), 6, i); 
    rowsChr = strcmp(T_out.Chromosome, uniqueChromosomes{i});
%     plotHiddenStatesOverTimeGW(T_out(rowsChr,:), daysLabels);
    plotHiddenStatesOverTime_TetKO(T_out(rowsChr,:), daysLabels);

    title(uniqueChromosomes{i});
    ylim([0,1]);
%     ylabel('Level');
%     legend('mm', 'toth', 'um-mu', 'uu');

end
print('-bestfit', '~/Desktop/levelPerChr_tet_a', '-dpdf');

fig = figure('Name', 'Levels per Chr',...
    'units', 'centimeters', ...
    'Position', [0, 0, 30, 30]);

for i=1:numOfChromosomes/2
    subplot(ceil(numOfChromosomes/12), 6, i); 
    rowsChr = strcmp(T_out.Chromosome, uniqueChromosomes{i+numOfChromosomes/2});
%     plotHiddenStatesOverTimeGW(T_out(rowsChr,:), daysLabels);
    plotHiddenStatesOverTime_TetKO(T_out(rowsChr,:), daysLabels);
    
    title(uniqueChromosomes{i+numOfChromosomes/2});
    ylim([0,1]);
%     ylabel('Level');
%     legend('mm', 'toth', 'um-mu', 'uu');

end
print('-bestfit', '~/Desktop/levelPerChr_tet_b', '-dpdf');
close all;


end