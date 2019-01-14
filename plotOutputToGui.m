%function [figure1] = plotOutput(x, p0, days, obsBis, obsOx, EBis, EOx, Cov, dataPointsName, dataLabels, process)
function plotOutputToGui(axesArray, xminm, pData, pModel, pAllStates,...
    CovPlot, dataPointsName, dataPoints, dataLabels, paramNames, knots, maxDegree, ubR)


%TO BE GENERALIZED FOR ANOTHER CASE
days = dataPoints;
numOfDays = size(days, 1);
numOfExperiments = size(pData, 3);

switch numOfExperiments
    case 1
        %define states' subsets
        uu = pAllStates(:,1);
        %um+mu
        um = sum([pAllStates(:,2), pAllStates(:,5)], 2);
        mm = pAllStates(:,6);
    case 2 
        %define states' subsets
        uu = pAllStates(:,1);
        %um+mu
        um = sum([pAllStates(:,2), pAllStates(:,5)], 2);
        %uh+hu
        uh = sum([pAllStates(:,3), pAllStates(:,9)], 2);
        %mh+hm
        hm = sum([pAllStates(:,7), pAllStates(:,10)], 2);
        hh = pAllStates(:,11);
        toth = hm + uh + hh;
        mm = pAllStates(:,6);
        hydroxyStates = [uh, hm, hh];
    case 3
        %define states' subsets
        uu = pAllStates(:,1);
        %um+mu
        um = sum([pAllStates(:,2), pAllStates(:,5)], 2);
        %uh+hu
        uh = sum([pAllStates(:,3), pAllStates(:,9)], 2);
        %mh+hm
        hm = sum([pAllStates(:,7), pAllStates(:,10)], 2);
        hh = pAllStates(:,11);
        toth = hm + uh + hh;
        mm = pAllStates(:,6);
        hydroxyStates = [uh, hm, hh];
        %uf+fu
        uf = sum([pAllStates(:,4), pAllStates(:,13)], 2);
        %mf+fm
        mf = sum([pAllStates(:,8), pAllStates(:,14)], 2);
        %hf+fh
        hf = sum([pAllStates(:,12), pAllStates(:,15)], 2);
        ff = pAllStates(:,16);
        totf = uf+mf+hf+ff;
        formalStates = [uf, mf, hf, ff];
end



%% --create axes1
% h = axes('Parent', fig, 'units', 'normalized', 'XTick', days',...
%     'Position', [0.31250000000000006 0.6045197740112994 0.1993243243243244 0.288135593220339]); 
h = axesArray(1);


h.XTick = days';
h.XLim =  [days(1), days(numOfDays)+0.1];
h.YLim =  [0, 1]; 
%fontSize as a proportion of the axes height (FontSize in pts is 14)
h.FontUnits = 'normalized';
h.FontSize = 0.0664;
set(h, 'FontUnits', 'Points');
fontSizePts = ceil(get(h, 'FontSize'));
%h.FontSize = 14;

%define colors you want to use
rgbcolor1 = [0         0.4470    0.7410;       
            rgb('PaleGreen');
            rgb('DarkGreen');
            0.8500    0.3250    0.0980];

% Uncomment the following line to preserve the Y-limits of the axes
ylim(h, [0 1]);
xlabel(h, dataPointsName, 'FontSize', floor(fontSizePts*1.3));
ylabel(h, 'frequency', 'FontSize', floor(fontSizePts*1.45));

h.TitleFontSizeMultiplier = 1.5;
% h.LabelFontSizeMultiplier = 1.5; %problem is this changes both xlabel and ylabel size
title(h, 'BS')

box(h, 'on');      
hold(h, 'on');

for i=1:4
    plot(h, days, pData(:,i, 1), 'Color', rgbcolor1(i,:), 'LineWidth', 2);
    hold(h, 'on');
    plot(h, days, pModel(:,i,1), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
    hold(h, 'on');
end

% box(h, 'on');        
% hold(h, 'on');

% axes(h);

labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
legend1 = legend(h, labels, 'Location', 'Best', 'FontSize', floor(fontSizePts/1.7)); 
set(legend1, 'ButtonDownFcn', [], 'FontUnits', 'normalized');



%% -- Create axes2
% h = axes('Parent', fig, 'Units', 'normalized', ...
%     'Position', [0.5363175675675675 0.6045197740112994 0.19847972972972971 0.28954802259887014]); 

h = axesArray(2);

%'XTickLabel', strsplit(int2str(days')), 'XTick', days', ...
%     'XLim', [days(1),  days(numOfDays)], 'YLim', [0, 1], 'FontSize', 12, ...
%     'Position', [0.5363175675675675 0.6045197740112994 0.19847972972972971 0.28954802259887014]);

% h.XTickLabel =  strsplit(int2str(days'));
h.XTick = days';
h.XLim =  [days(1), days(numOfDays)+0.1];
h.YLim =  [0, 1];

%   h.FontSize = 14;
%fontSize as a proportion of the axes height (FontSize in pts is 14)
h.FontUnits = 'normalized';
h.FontSize = 0.0664;
set(h, 'FontUnits', 'Points');
fontSizePts = ceil(get(h, 'FontSize'));

% Uncomment the following line to preserve the Y-limits of the axes
ylim(h, [0 1]);   
h.TitleFontSizeMultiplier = 1.5;
title(h, 'oxBS')
%   title(h, 'oxBS', 'FontSize', 20)
xlabel(h, dataPointsName, 'FontSize', floor(fontSizePts*1.3));
%   xlabel(h, handles.dataPointsName, 'FontSize', 18);
box(h, 'on');        
hold(h, 'on');

%define colors you want to use
rgbcolor1 = [0         0.4470    0.7410;       
            rgb('PaleGreen');
            rgb('DarkGreen');
            0.8500    0.3250    0.0980];

                
if numOfExperiments > 1
    
    for i=1:4
        plot(h, days, pData(:,i,2), 'Color', rgbcolor1(i,:), 'LineWidth',2);
        hold(h, 'on');
        plot(h, days, pModel(:,i,2), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
        hold(h, 'on');
    end
    
    labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
    legend1 = legend(h, labels, 'Location', 'Best', 'FontSize', floor(fontSizePts/1.7)); 
    set(legend1, 'ButtonDownFcn', [], 'FontUnits', 'normalized', 'Visible', 'off');   

    

else
    %set(dHandles.annot1, 'Visible', 'on');
end    


%% -- Create axes3
% h = axes('Parent', fig, 'Units', 'normalized', ...
%     'Position', [0.5363175675675675 0.6045197740112994 0.19847972972972971 0.28954802259887014]); 

h = axesArray(3);

%'XTickLabel', strsplit(int2str(days')), 'XTick', days', ...
%     'XLim', [days(1),  days(numOfDays)], 'YLim', [0, 1], 'FontSize', 12, ...
%     'Position', [0.5363175675675675 0.6045197740112994 0.19847972972972971 0.28954802259887014]);

% h.XTickLabel =  strsplit(int2str(days'));
h.XTick = days';
h.XLim =  [days(1), days(numOfDays)+0.1];
h.YLim =  [0, 1];

%   h.FontSize = 14;
%fontSize as a proportion of the axes height (FontSize in pts is 14)
h.FontUnits = 'normalized';
h.FontSize = 0.0664;
set(h, 'FontUnits', 'Points');
fontSizePts = ceil(get(h, 'FontSize'));

% Uncomment the following line to preserve the Y-limits of the axes
ylim(h, [0 1]);   
h.TitleFontSizeMultiplier = 1.5;
title(h, 'mabBS')
%   title(h, 'oxBS', 'FontSize', 20)
xlabel(h, dataPointsName, 'FontSize', floor(fontSizePts*1.3));
%   xlabel(h, handles.dataPointsName, 'FontSize', 18);
box(h, 'on');        
hold(h, 'on');

%define colors you want to use
rgbcolor1 = [0         0.4470    0.7410;       
            rgb('PaleGreen');
            rgb('DarkGreen');
            0.8500    0.3250    0.0980];
               
if (numOfExperiments > 2)
    
    for i=1:4
        plot(h, days, pData(:,i,3), 'Color', rgbcolor1(i,:), 'LineWidth',2);
        hold(h, 'on');
        plot(h, days, pModel(:,i,3), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
        hold(h, 'on');
    end
    
    labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
    legend1 = legend(h, labels, 'Location', 'Best', 'FontSize', floor(fontSizePts/1.7)); 
    set(legend1, 'ButtonDownFcn', [], 'FontUnits', 'normalized', 'Visible', 'off');   

    

else
    %set(dHandles.annot1, 'Visible', 'on');
end    






%% Create axes4
% h = axes('Parent',fig, 'Units', 'normalized',...
%     'Position', [0.31250000000000006 0.1807909604519774 0.42314189189189194 0.3545197740112994]);
h = axesArray(4);

%create xTickLabel
h.XTickLabel = dataLabels;
h.XTick = 1:numOfDays;

%fontSize as a proportion of the axes height (FontSize in pts is 12)
h.FontUnits = 'normalized';
h.FontSize = 0.0564;
set(h, 'FontUnits', 'Points');
fontSizePts = ceil(get(h, 'FontSize'));

% Uncomment the following line to preserve the Y-limits of the axes
ylim(h, [0 1]);
xlim(h, [0, numOfDays+1]);
box(h,'on');
hold(h,'on');

% Create ylabel
ylabel(h, 'level per states', 'HorizontalAlignment', 'center', 'FontSize', floor(fontSizePts*1.45));

if numOfExperiments == 3
    plotedStates = [mm, toth, um, totf, uu];
    % Create multiple lines using matrix input to bar
    bar1 = bar(plotedStates, 'Parent', h, 'BarLayout','stacked');
    set(bar1(1),'DisplayName','fullymethylated','FaceColor', [0.85 0.325 0.098]);
    set(bar1(2),'DisplayName','hydroxylated','FaceColor', [0.929 0.694 0.125]);
    set(bar1(3),'DisplayName','hemimethylated','FaceColor', [0.466 0.674 0.188]);
    set(bar1(4),'DisplayName','formal','FaceColor', rgb('LightSkyBlue'));
    set(bar1(5),'DisplayName','unmethylated','FaceColor', [0 0.447 0.741]);
elseif numOfExperiments == 2
    plotedStates = [mm, toth, um, uu];
    % Create multiple lines using matrix input to bar
    bar1 = bar(plotedStates, 'Parent', h, 'BarLayout','stacked');
    set(bar1(1),'DisplayName','fullymethylated','FaceColor', [0.85 0.325 0.098]);
    set(bar1(2),'DisplayName','hydroxylated','FaceColor', [0.929 0.694 0.125]);
    set(bar1(3),'DisplayName','hemimethylated','FaceColor', [0.466 0.674 0.188]);
    set(bar1(4),'DisplayName','unmethylated','FaceColor', [0 0.447 0.741]);
    
else
    % Create multiple lines using matrix input to bar
    plotedStates = [mm, um, uu];
    bar1 = bar(plotedStates, 'Parent', h, 'BarLayout','stacked');
    set(bar1(1),'DisplayName','fullymethylated','FaceColor', [0.85 0.325 0.098]);
    set(bar1(2),'DisplayName','hemimethylated','FaceColor', [0.466 0.674 0.188]);
    set(bar1(3),'DisplayName','unmethylated','FaceColor', [0 0.447 0.741]);
    
end

% Create legend
legend4 = legend(h, 'Location', 'best');
set(legend4, 'ButtonDownFcn',[], 'FontSize', ceil(fontSizePts*0.7));    





%% Create axes5
% h = axesArray('Parent',fig, 'Units', 'normalized',...
%     'Position',[0.7761824324324325 0.1807909604519774 0.21875 0.3545197740112994]);

h = axesArray(5);

h.XTickLabel = dataLabels;
h.XTick = 1:numOfDays;

h.FontUnits = 'normalized';
h.FontSize = 0.04;
set(h, 'FontUnits', 'Points');
fontSizePts = ceil(get(h, 'FontSize'));

% h.FontSize = 10;

colormap(h, [0.929 0.5 0.125;0.929 0.694 0.125;0.929 0.9 0.125]);
%Uncomment the following line to preserve the X-limits of the axes
h.XLim = ([0 numOfDays+1]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(h,[0, max(sum(hydroxyStates, 2)) + 0.05]);
h.YLim = ([0, 0.5]); %as in the paper
box(h,'on');
hold(h,'on');

% Create ylabel
ylabel(h, 'hydroxylation level', 'HorizontalAlignment', 'center', 'FontSize', ceil(fontSizePts*1.6));

if numOfExperiments > 1
    % Create multiple lines using matrix input to bar
    bar2 = bar(hydroxyStates,'Parent', h, 'BarLayout', 'stacked');
    set(bar2(1),'DisplayName','hm-mh', 'FaceColor', [0.929 0.5 0.125]);
    set(bar2(2),'DisplayName','uh-hu', 'FaceColor', [0.929 0.694 0.125]);
    set(bar2(3),'DisplayName','hh', 'FaceColor', rgb('Yellow'));

    % Create legend
    legend5 = legend(h, 'Location', 'northwest'); 
    set(legend5, 'ButtonDownFcn',[]);
end


%% Create axes6
% h = axesArray('Parent',fig, 'Units', 'normalized',...
%     'Position',[0.7761824324324325 0.1807909604519774 0.21875 0.3545197740112994]);

h = axesArray(6);

h.XTickLabel = dataLabels;
h.XTick = 1:numOfDays;

h.FontUnits = 'normalized';
h.FontSize = 0.04;
set(h, 'FontUnits', 'Points');
fontSizePts = ceil(get(h, 'FontSize'));

% h.FontSize = 10;

% colormap(h, [0.929 0.5 0.125; 0.929 0.694 0.125; 0.929 0.9 0.125]);
%Uncomment the following line to preserve the X-limits of the axes
h.XLim = ([0 numOfDays+1]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(h,[0, max(sum(formalStates, 2)) + 0.1]);
h.YLim = ([0, 0.5]); %as in the paper
box(h,'on');
hold(h,'on');

% Create ylabel
ylabel(h, '5fC level', 'HorizontalAlignment', 'center', 'FontSize', ceil(fontSizePts*1.6));

if numOfExperiments > 2
    % Create multiple lines using matrix input to bar
    bar3 = bar(formalStates,'Parent', h, 'BarLayout', 'stacked');
    set(bar3(1),'DisplayName','mf-fm', 'FaceColor', rgb('LightSkyBlue'));
    set(bar3(2),'DisplayName','uf-fu', 'FaceColor', rgb('DeepSkyBlue'));
    set(bar3(3),'DisplayName','hf-fh', 'FaceColor', rgb('DodgerBlue'));
    set(bar3(4),'DisplayName','ff', 'FaceColor', rgb('CornflowerBlue'));

    % Create legend
    legend6 = legend(h, 'Location', 'northwest'); 
    set(legend6, 'ButtonDownFcn',[]);
end

%% Create axes7
% h = axes('Parent', fig, 'Units', 'normalized', ...
%     'Position',[0.7753378378378377 0.6031073446327684 0.22043918918918903 0.2853107344632768]);

h = axesArray(7);


box(h, 'on');
hold(h, 'on');


t = days;
%names to variables
x = xminm;
Cov = CovPlot;

% h.FontSize = 14;
%fontSize as a proportion of the axes height (FontSize in pts is 14)
h.FontUnits = 'normalized';
h.FontSize = 0.0664;
set(h, 'FontUnits', 'Points');
fontSizePts = ceil(get(h, 'FontSize'));

% Create ylabel
ylabel(h, 'efficiency', 'FontSize', floor(fontSizePts*1.45));
% Create xlabel
xlabel(h, dataPointsName, 'FontSize', floor(fontSizePts*1.3));


%the next is to be put as an argument above if we want specific Ticks
xlim(h, [t(1)-0.08, t(numOfDays)+0.08]);

% ylim('auto');
% y = ylim;
% ylim(h, [0, y(2)]);
box(h, 'on');
hold(h, 'on');

if ~isempty(x)

    %hydroxylation probability is included or not
    if strcmp(paramNames{end}, 'p')
        numOfRateParams = size(paramNames, 2)-1;
    else
        numOfRateParams = size(paramNames, 2);
    end    

    %find the degree (or numOfCoef) of the underlying polynom
    % numOfCoef = ind(2) - 1;
    numOfCoef = maxDegree + 1;
    numOfRates = floor(numOfRateParams/numOfCoef);
    ratesNames = cell(1, numOfRates);
    %get the rates names from paramNames vector
    for i=1:numOfRates
        ratesNames{i} = paramNames{(i-1)*numOfCoef+1};
    end
    
    allRatesArray = {'maint', 'deNovo', 'hydr', 'form', 'act_dem'};
    allLabelsArray = {'\mu_m', '\mu_d', '\eta', '\phi', '\delta'};

    [~, idb] = ismember(ratesNames, allRatesArray);
    
    labels = allLabelsArray(nonzeros(idb));

    rgbColor = [0.8500    0.3250    0.0980; %Red
                 0         0.4470    0.7410;  %DarkBlue
                 0.9290    0.6940    0.1250;  %yellow
                 0         0.20    0.5410;      %lightblue
                 rgb('Teal');                   %green
                 0.5430    0         0   ];     %DarkRed

    %define and plot each rate as a function of time
    numOfPlottedRates = size(labels, 2);
    
    tPlot = 0:0.5:days(end);
    l = zeros(1, numOfPlottedRates);
    p = zeros(1, numOfPlottedRates);

    [y, varOfy] = getEfficiencies(x, maxDegree, tPlot, knots, Cov);
    y = y(:, nonzeros(idb));
    varOfy = varOfy(:, nonzeros(idb));
    
    axes(h);
    for i=1:numOfPlottedRates
        [l(i),p(i)] = boundedline(tPlot, y(:,i), 1.96*sqrt(varOfy(:,i)), 'alpha', 'cmap', rgbColor(idb(i),:));
        %     outlinebounds(l,p);
        hold(h, 'on');
    end
    
    ylim(h, [0, ubR]);

    legend(l, labels, 'Location', 'east', 'FontSize', 12);
    legend('boxoff')

end    
% legend7 = legend(h, l, labels, 'Location', 'east', 'FontSize', ceil(fontSizePts*0.85));
% % legend7.String = legend7.String(2:2:size(legend7.String));
% set(legend7, 'ButtonDownFcn', []);

%the next is to be put as argument above if we want specific Ticks
% xlim(h, [t(1)-0.08, t(numOfDays)+0.08]);
% ylim(h, [0, 1]);
% box(h, 'on');
% hold(h, 'on');

% for i=1:2:2*numOfRates-1
%     theta((i+1)/2,:) = x(i) + x(i+1)*t;
%     %b0 should equal b1 for every parameter
%     %theta((i+1)/2,1) = x(i) + x(i+1)*1;
%     var_b0 = Cov(i,i);
%     var_b1 = Cov(i+1,i+1);
%     cov = Cov(i,i+1);
%     var_theta((i+1)/2,:) = var_b0 + var_b1*t.^2 + 2*cov*t;
%     errorbar(t, theta((i+1)/2,:), sqrt(var_theta((i+1)/2,:)), 'Parent', h, 'Marker', 'o', 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2((i+1)/2,:))
%     hold(h, 'on');
% end
% %store and plot the lambda = mu1+mu2-mu1*mu2 over time
% %theta(k,:) is the k-th entry of v =[mu_m, mu_d, eta] 
% covMu_mMu_d(1,:) = Cov(1,3) + t.*Cov(1,4) + t.*Cov(2,3) + t.^2.*Cov(2,4);
% %some necessary moments
% % secondMu_m(1,:) = var_theta(1,:) + theta(1,:).^2;
% % secondMu_d(1,:) = var_theta(2,:) + theta(2,:).^2;
% secondMu_mMu_d(1,:) = covMu_mMu_d(1,:) + theta(1,:).*theta(2,:);
% thirdMu_m2Mu_d(1,:) = theta(1,:).^2.*theta(2,:) + 2*covMu_mMu_d(1,:).*theta(1,:) + var_theta(1,:).*theta(2,:);
% thirdMu_mMu_d2(1,:) = theta(2,:).^2.*theta(1,:) + 2*covMu_mMu_d(1,:).*theta(2,:) + var_theta(2,:).*theta(1,:);
% fourthMu_m2Mu_d2(1,:) = theta(1,:).^2 .*theta(2,:).^2 + var_theta(1,:).*var_theta(2,:) + 2*covMu_mMu_d(1,:).^2 + ...
%     4*covMu_mMu_d(1,:).*theta(1,:).*theta(2,:) + var_theta(1,:).*theta(2,:).^2 + var_theta(2,:).*theta(1,:).^2;
% 
% %the total methylation of hemimethylated sites (lambda)
% theta(numOfRates+1,:) = theta(1,:) + theta(2,:) - theta(1,:).*theta(2,:);
% 
% %quantities refer to the appendix
% varMu_mMu_d(1,:) = fourthMu_m2Mu_d2(1,:) - secondMu_mMu_d(1,:).^2;
% CovMu_m2Mu_d(1,:) = thirdMu_m2Mu_d(1,:) - theta(1,:).*secondMu_mMu_d(1,:);
% CovMu_mMu_d2(1,:) = thirdMu_mMu_d2(1,:) - theta(2,:).*secondMu_mMu_d(1,:);
% 
% %sum all the previous to get the variance of lambda
% var_theta(numOfRates+1,:) = var_theta(1,:) + var_theta(2,:) + varMu_mMu_d(1,:) + 2*(covMu_mMu_d(1,:) - CovMu_m2Mu_d(1,:) - CovMu_mMu_d2(1,:));
% 
% %plot the standard deviation of lambda
% errorbar(t, theta(numOfRates+1,:), sqrt(var_theta(numOfRates+1,:)), 'Parent', h, 'Marker', 'o', 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2(numOfRates+1,:));

% legend7 = legend(h, l, labels, 'Location', 'east', 'FontSize', ceil(fontSizePts*0.85));
% % legend7.String = legend7.String(2:2:size(legend7.String));
% set(legend7, 'ButtonDownFcn', []);



end
