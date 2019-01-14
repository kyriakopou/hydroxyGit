%function [figure1] = plotOutput(x, p0, days, obsBis, obsOx, EBis, EOx, Cov, dataPointsName, dataLabels, process)
function [figure1] = plotOutput(x, pData, pModel, pAllStates, Cov, dataPointsName, days, dataLabels, paramNames, knots, maxDegree, ubR)


%TO BE GENERALIZED FOR ANOTHER CASE
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


%call DSHydroxy with the optimum parameter values
%[~, ~, pBis, pOx, pBisModel, pOxModel, pAllStates] = DSHydroxyEmbryo(x, p0, days, obsBis, obsOx, EBis, EOx, process);

%display the predicted level of all states over time
% fprintf('Predicted level of all states \n')
% fprintf('\tuu \t um \t mu \t\t uh  \t  hu  \t  hm  \t  mh  \t  mm  \t  hh\n');
% display(pAllStates);


%display and plot the normalized level  taken from the data of all states over time
%using MLEstimator
% pData2 = zeros(numOfDays, 9);
% for t=1:numOfDays
%     pData2(t,:) = maxLikelihood(obsBis(t,:), obsOx(t,:), EBis(:,:,t), EOx(:,:,t)); 
% end
% 
% % fprintf('\n Level of all states taken from data by MLEstimator \n')
% % fprintf('\tuu \t um \t mu \t\t uh  \t  hu  \t  hm  \t  mh  \t  mm  \t  hh\n');
% % display(pData2)
% uu = pData2(:,1);
% um = sum(pData2(:,2:3), 2);
% uh = sum(pData2(:,4:5), 2);
% hm = sum(pData2(:,6:7), 2);
% hh = pData2(:,9);
% toth = hm + uh + hh;
% mm = pData2(:,8);
% hydroxyStates = [uh, hm, hh];
% plotedStates = [mm, toth, um, uu];
% createFigurePaper(plotedStates, hydroxyStates, days, regionName);




%plot all states over time
% figure('name', 'Data vs Model');
% plot(days, pAllStates);
% labels = {'uu', 'um', 'mu', 'uh', 'hu', 'hm', 'mh', 'mm', 'hh' };
% xlabel('days', 'FontSize', 20);
% ylabel('frequency','FontSize',20);
% plegend = legend(labels, 'Location', 'BestOutside', 'Orientation', 'vertical');
% set(plegend, 'Fontsize', 14)




% Create figure and fit the window size in Position (IMPORTANT)
% figure1 = figure('visible', 'on', 'Name','Model prediction',...
%     'units', 'normalized', ...
%     'Colormap',[0.929 0.5 0.125;0.929 0.694 0.125;0.929 0.9 0.125], ...
%     'Position', [0, 0, 0.9, 0.75]);
% 
% set(figure1, 'Units', 'centimeters');
% fontSizePts = ceil(get(figure1, 'Position'));

figure1 = figure('visible', 'off', 'Name', 'Model prediction',...
    'units', 'centimeters', ...
    'Colormap',[0.929 0.5 0.125;0.929 0.694 0.125;0.929 0.9 0.125], ...
    'Position', [0, 0, 41, 22]);

    
%define colors you want to use
rgbcolor1 = [0         0.4470    0.7410;       
            rgb('PaleGreen');
            rgb('DarkGreen');
            0.8500    0.3250    0.0980];


%% Create axes1
axes1 = axes('Parent', figure1, 'XTickLabel', strsplit(int2str(days')), 'XTick', days', ...
    'XLim', [days(1),  days(numOfDays)], 'YLim', [0, 1], 'FontSize', 12, ...
    'Position', [0.143767361111111 0.59006379585327 0.167951388888889 0.281180223285486]);

box(axes1,'on');
hold(axes1,'on');            
        
title('BS', 'FontSize', 15)

xlabel(dataPointsName, 'FontSize', 15);
% yLabelName = strcat('frequency (', regionName,')');
yLabelName = 'frequency';
ylabel(yLabelName, 'FontSize', 15);
    
     
for i=1:4
    plot(days, pData(:,i,1), 'Color', rgbcolor1(i,:), 'LineWidth',2);
    hold on;
    plot(days, pModel(:,i,1), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
    hold on;
end

% labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
% leg1 = legend(labels, 'location', 'northeast', 'fontsize', 8);

labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
% set(figOut, 'CurrentAxes', h(1));
% columnlegend(2, labels, 'fontsize', 7);
legendflex(labels, 'anchor', {'nw','nw'}, ...
        'buffer', [-5 -5], ...
        'ncol', 4, ...
        'fontsize', 6, ...
        'xscale', 0.6, ...
        'box', 'off');


%% Create axes2
axes2 = axes('Parent', figure1, 'XTickLabel', strsplit(int2str(days')), 'XTick', days', ...
    'XLim', [days(1),  days(numOfDays)], 'YLim', [0, 1], 'FontSize', 12, ...
    'Position', [0.3340 0.59006379585327 0.16951388888889 0.281180223285486]);


box(axes2,'on');
hold(axes2,'on');

title('oxBS', 'FontSize', 15)
xlabel(dataPointsName, 'FontSize', 15);

if numOfExperiments >= 2

    for i=1:4
        plot(days, pData(:,i,2), 'Color', rgbcolor1(i,:), 'LineWidth', 2);
        hold on;
        plot(days, pModel(:,i,2), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
        hold on;
    end
    
%     labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
%     leg2 = legend(labels, 'location', 'northeast', 'fontsize', 8);
       
end


%% Create axes3
axes3 = axes('Parent', figure1, 'XTickLabel', strsplit(int2str(days')), 'XTick', days', ...
    'XLim', [days(1),  days(numOfDays)], 'YLim', [0, 1], 'FontSize', 12, ...
    'Position', [0.5280 0.59006379585327 0.16951388888889 0.281180223285486]);


box(axes3,'on');
hold(axes3,'on');

title('mabBS', 'FontSize', 15)
xlabel(dataPointsName, 'FontSize', 15);

if numOfExperiments == 3

    for i=1:4
        plot(days, pData(:,i,3), 'Color', rgbcolor1(i,:), 'LineWidth', 2);
        hold on;
        plot(days, pModel(:,i,3), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
        hold on;
    end
    
%     labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
%     leg2 = legend(labels, 'location', 'northeast', 'fontsize', 8);
       
end



%% Create axes
axes4 = axes('Parent', figure1, 'FontSize', 10,...
    'XTick', 0:days(numOfDays), ...
    'Position',[0.74084375 0.59006379585327 0.169582849399922 0.281180223285485]);

box(axes4,'on');
hold(axes4,'on');

if ~isempty(x)
    %variables names
    t = days;
    numOfDays = size(t,1);

    % Create ylabel
    % yLabelName = strcat('efficiency (', regionName,')');
    yLabelName = 'efficiency';
    ylabel(yLabelName, 'FontSize', 14);
    % Create xlabel
    xlabel(dataPointsName, 'FontSize', 15);

    %hydroxylation probability is included or not
    if strcmp(paramNames{end}, 'p')
        numOfParams = size(x, 2)-1;
    else
        numOfParams = size(x, 2);
    end    

    %find the degree (or numOfCoef) of the underlying polynom
    numOfCoef = maxDegree + 1;

    numOfRates = floor(numOfParams / numOfCoef);
    ratesNames = cell(1, numOfRates);
    %get the rates names from paramNames vector
    for i=1:numOfCoef:numOfParams
        ratesNames{floor(i/numOfCoef) + 1} = paramNames{i}(4:end);
    end

    allLabelsArray = {'\mu_m', '\mu_d', '\eta', '\phi', '\delta'};
    % [~, idb] = ismember(ratesNames, allRatesArray);
    % labels = allLabelsArray(idb);

    rgbColor = [0.8500    0.3250    0.0980; %Red
                 0         0.4470    0.7410;  %DarkBlue
                 0.9290    0.6940    0.1250;  %yellow
                 0         0.20    0.5410;      %lightblue
                 rgb('Teal');                   %green
                 0.5430    0         0   ];     %DarkRed

    %compute the degrees of freedom (dof)

    %define and plot each rate as a function of time
    theta = zeros(numOfRates+1, numOfDays);
    var_theta = zeros(numOfRates+1, numOfDays);


    %the next is to be put as argument above if we want specific Ticks
    xlim(axes4, [t(1)-0.08, t(numOfDays)+0.08]);
    ylim(axes4, [0, ubR]);
    box(axes4, 'on');
    hold(axes4, 'on');


    tPlot = 0:0.5:days(end);

    l = zeros(1, numOfRates);
    p = zeros(1, numOfRates);


    [y, varOfy] = getEfficiencies(x, maxDegree, tPlot, knots, Cov);

%   find(strcmp(ratesNames(i), allRatesArray))
    for i=1:numOfRates
        [l(i), p(i)] = boundedline(tPlot, y(:,i), 1.96*sqrt(varOfy(:,i)), 'alpha', 'cmap', rgbColor(i,:));
        %     outlinebounds(l,p);
        hold(axes4, 'on');
    end
    
    ylim(axes4, [0, ubR]);

    legend(l, allLabelsArray, 'Location', 'east', 'FontSize', 12);
    legend('boxoff');
    
end
% maxDegree = numOfCoef -1;
% [~, mu2, h, f, dem] = getEfficiencies(x, maxDegree, tPlot, knots);


% for i=1:numOfCoef:numOfRateParams-1
%     theta((i+1)/numOfCoef,:) = x(i) + x(i+1)*t;
%     %b0 should equal b1 for every parameter
%     %theta((i+1)/2,1) = x(i) + x(i+1)*1;
%     var_b0 = Cov(i,i);
%     var_b1 = Cov(i+1,i+1);
%     cov = Cov(i,i+1);
%     var_theta((i+1)/2,:) = var_b0 + var_b1*t.^2 + 2*cov*t;
%     errorbar(t, theta((i+1)/2,:), sqrt(var_theta((i+1)/2,:)), 'Marker', 'o', 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2((i+1)/2,:))
%     hold(axes4, 'on');
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
% errorbar(t, theta(numOfRates+1,:), sqrt(var_theta(numOfRates+1,:)), 'Marker', 'o', 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2(numOfRates+1,:));


% set(legend5,...
%     'Position', [0.654281255142446 0.413078149920255 0.0453125 0.106858054226475],...   
%     'FontSize', 12);


%% Create axes5
axes5 = axes('Parent', figure1, 'XTickLabel', dataLabels, ...
    'XTick', 1:numOfDays,...
    'FontSize', 12,...
    'Position', [0.143767361111111 0.102840909090909 0.27875 0.400715709728868]);
%'Position',[0.206736353077816 0.102840909090909 0.372822299651568 0.7]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes5,[0 1]);
box(axes5,'on');
hold(axes5,'on');


% Create ylabel
% yLabelName = strcat('level per states (', regionName, ')');
yLabelName = 'level per states ';
ylabel(yLabelName, 'FontSize', 14);
if numOfExperiments == 3
    plotedStates = [mm, toth, um, totf, uu];
    % Create multiple lines using matrix input to bar
    bar1 = bar(plotedStates, 'BarLayout','stacked');
    set(bar1(1),'DisplayName','fullymethylated','FaceColor', [0.85 0.325 0.098]);
    set(bar1(2),'DisplayName','hydroxylated','FaceColor', [0.929 0.694 0.125]);
    set(bar1(3),'DisplayName','hemimethylated','FaceColor', [0.466 0.674 0.188]);
    set(bar1(4),'DisplayName','formal','FaceColor', rgb('LightSkyBlue'));
    set(bar1(5),'DisplayName','unmethylated','FaceColor', [0 0.447 0.741]);
elseif numOfExperiments == 2
    plotedStates = [mm, toth, um, uu];
    % Create multiple lines using matrix input to bar
    bar1 = bar(plotedStates, 'BarLayout','stacked');
    set(bar1(1),'DisplayName','fullymethylated','FaceColor', [0.85 0.325 0.098]);
    set(bar1(2),'DisplayName','hydroxylated','FaceColor', [0.929 0.694 0.125]);
    set(bar1(3),'DisplayName','hemimethylated','FaceColor', [0.466 0.674 0.188]);
    set(bar1(4),'DisplayName','unmethylated','FaceColor', [0 0.447 0.741]);
    
else
    % Create multiple lines using matrix input to bar
    plotedStates = [mm, um, uu];
    bar1 = bar(plotedStates, 'BarLayout', 'stacked');
    set(bar1(1),'DisplayName','fullymethylated','FaceColor', [0.85 0.325 0.098]);
    set(bar1(2),'DisplayName','hemimethylated','FaceColor', [0.466 0.674 0.188]);
    set(bar1(3),'DisplayName','unmethylated','FaceColor', [0 0.447 0.741]);
    
end

% Create legend
% leg4 = legend(axes4, 'Location', 'best');
% set(leg4, 'FontSize', 10);


%% Create axes6
axes6 = axes('Parent', figure1, 'XTickLabel', dataLabels,...
    'XTick',1:numOfDays,...
    'FontSize', 10,...
    'Position',[0.46484375 0.102840909090909 0.199582849399922 0.400715709728868]);


%Uncomment the following line to preserve the X-limits of the axes
xlim(axes6,[0 numOfDays+1]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes6,[0, 0.50]); %as in the paper
ylim(axes6,[0, max(sum(hydroxyStates, 2)) + 0.05]);
box(axes6,'on');
hold(axes6,'on');

% Create ylabel
ylabel('hydroxylation level', 'FontSize', 14);

if numOfExperiments > 1
    % Create multiple lines using matrix input to bar
    bar2 = bar(hydroxyStates,'Parent', axes6, 'BarLayout', 'stacked');
    set(bar2(1),'DisplayName', 'hm-mh', 'FaceColor',  [0.929 0.5 0.125]);
    set(bar2(2),'DisplayName', 'uh-hu', 'FaceColor', [0.929 0.694 0.125]);
    set(bar2(3),'DisplayName', 'hh', 'FaceColor', rgb('Yellow'));
    
    % Create legend
    legend5 = legend(axes6, 'Location', 'northwest'); 
    set(legend5, 'ButtonDownFcn',[]);
end

%% Create axes7
axes7 = axes('Parent', figure1, 'XTickLabel', dataLabels,...
    'XTick',1:numOfDays,...
    'FontSize', 10,...
    'Position',[0.710484375 0.102840909090909 0.199582849399922 0.400715709728868]);


%Uncomment the following line to preserve the X-limits of the axes
xlim(axes7,[0 numOfDays+1]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes7,[0, 0.50]); %as in the paper
%ylim(axes5,[0, max(sum(hydroxyStates, 2)) + 0.05]);
box(axes7,'on');
hold(axes7,'on');

% Create ylabel
ylabel('5fC level', 'FontSize', 14);

if numOfExperiments == 3
    % Create multiple lines using matrix input to bar
    bar2 = bar(formalStates,'Parent', axes7, 'BarLayout', 'stacked');
    set(bar2(1),'DisplayName','mf-fm', 'FaceColor', rgb('LightSkyBlue'));
    set(bar2(2),'DisplayName','uf-fu', 'FaceColor', rgb('DeepSkyBlue'));
    set(bar2(3),'DisplayName','hf-fh', 'FaceColor', rgb('DodgerBlue'));
    set(bar2(4),'DisplayName','ff', 'FaceColor', rgb('CornflowerBlue'));

    % Create legend
    legend5 = legend(axes7, 'Location', 'northwest'); 
    set(legend5, 'ButtonDownFcn',[]);
end




end
