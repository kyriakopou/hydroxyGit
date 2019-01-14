function [figure1] = plotOutputSingleCpG(x, pBis, pOx, pBisModel, pOxModel, pAllStates, Cov, dataPointsName, days, dataLabels, dataType, regionName)
%IT SHOULD BE CALLED ONLY TO GET PDF OUTPUT FOR SINGLE CpG data (so that you get the 
%right distances and sizes of the plots)

numOfDays = size(days, 1);

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


uu = pAllStates(:,1);
um = sum(pAllStates(:,2:3), 2);
uh = sum(pAllStates(:,4:5), 2);
hm = sum(pAllStates(:,6:7), 2);
hh = pAllStates(:,9);
toth = hm + uh + hh;
mm = pAllStates(:,8);
hydroxyStates = [uh, hm, hh];
plotedStates = [mm, toth, um, uu];


% Create figure and fit the window size in Position (IMPORTANT)
figure1 = figure('visible', 'off', 'Name','Model prediction',...
    'units', 'normalized', ...
    'Colormap',[0.929 0.5 0.125;0.929 0.694 0.125;0.929 0.9 0.125], ...
    'Position', [0, 0, 0.9, 0.75]);

    
%define colors you want to use
rgbcolor1 = [0         0.4470    0.7410;       
            rgb('PaleGreen');
            rgb('DarkGreen');
            0.8500    0.3250    0.0980];


%% Create axes
axes1 = axes('Parent', figure1, 'XTickLabel', strsplit(int2str(days')), 'XTick', days', ...
    'XLim', [days(1),  days(numOfDays)], 'YLim', [0, 1], 'FontSize', 12, ...
    'Position', [0.183767361111111 0.59006379585327 0.167951388888889 0.281180223285486]);

box(axes1,'on');
hold(axes1,'on');            
        
title('BS', 'FontSize', 15)
xlabel('days', 'FontSize', 15);
yLabelName = strcat('frequency (', regionName,')');
ylabel(yLabelName, 'FontSize', 15);
    
     
for i=1:4
    plot(days, pBis(:,i), 'Color', rgbcolor1(i,:), 'LineWidth',2);
    hold on;
    plot(days, pBisModel(:,i), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
    hold on;
end

% labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
% leg1 = legend(labels, 'location', 'northeast', 'fontsize', 8);




%% Create axes
axes2 = axes('Parent', figure1, 'XTickLabel', strsplit(int2str(days')), 'XTick', days', ...
    'XLim', [days(1),  days(numOfDays)], 'YLim', [0, 1], 'FontSize', 12, ...
    'Position', [0.3740 0.59006379585327 0.16951388888889 0.281180223285486]);


box(axes2,'on');
hold(axes2,'on');

title('oxBS', 'FontSize', 15)
xlabel('days', 'FontSize', 15);

if (strcmp(dataType, 'oxBS'))

    for i=1:4
        plot(days, pOx(:,i), 'Color', rgbcolor1(i,:), 'LineWidth', 2);
        hold on;
        plot(days, pOxModel(:,i), 'Color', rgbcolor1(i,:), 'Marker' , 'o', 'LineStyle', '--', 'LineWidth', 2);
        hold on;
    end
    
%     labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
%     leg2 = legend(labels, 'location', 'northeast', 'fontsize', 8);
       
end



%% Create axes
axes3 = axes('Parent', figure1, 'FontSize', 12,...
    'XTick', 0:days(numOfDays), ...
    'Position',[0.61484375 0.59006379585327 0.189582849399922 0.281180223285485]);

box(axes3,'on');
hold(axes3,'on');

%variables names
t = days;
numOfDays = size(t,1);

% Create ylabel
if strcmp(regionName(end), int2str(1))
    yLabelName = 'efficiency';
else
    yLabelName = '';
end

ylabel(yLabelName, 'FontSize', 16);
% Create xlabel
xlabel('days', 'FontSize', 16);


%hydroxylation rate is included or not
if (strcmp(dataType, 'oxBS'))
    numOfRates = 3;
    labels = {'\mu_m', '\mu_d', '\eta', '\lambda'};
    rgbcolor2 = [0.8500    0.3250    0.0980;
             0         0.4470    0.7410;  %Red
             0.9290    0.6940    0.1250;  %yellow
             0.5430    0         0     ]; %DarkRed
else    
    numOfRates = 2;
    labels = {'\mu_m', '\mu_d', '\lambda'};
    rgbcolor2 = [0.8500    0.3250    0.0980;
             0         0.4470    0.7410;
             0.5430    0         0    ];  %DarkRed
end

%define and plot each rate as a function of time
theta = zeros(numOfRates+1, numOfDays);
var_theta = zeros(numOfRates+1, numOfDays);


% %the next is to be put as argument above if we want specific Ticks
%
xlim(axes3, [t(1)-0.08, t(numOfDays)+0.08]);
ylim(axes3, [0, 1]);
box(axes3, 'on');
hold(axes3, 'on');

for i=1:2:2*numOfRates-1
    theta((i+1)/2,:) = x(i) + x(i+1)*t;
    %b0 should equal b1 for every parameter
    %theta((i+1)/2,1) = x(i) + x(i+1)*1;
    var_b0 = Cov(i,i);
    var_b1 = Cov(i+1,i+1);
    cov = Cov(i,i+1);
    var_theta((i+1)/2,:) = var_b0 + var_b1*t.^2 + 2*cov*t;
    errorbar(t, theta((i+1)/2,:), sqrt(var_theta((i+1)/2,:)), 'Marker', 'o', 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2((i+1)/2,:))
    hold on;
end
%store and plot the lambda = mu1+mu2-mu1*mu2 over time
%theta(k,:) is the k-th entry of v =[mu_m, mu_d, eta] 
covMu_mMu_d(1,:) = Cov(1,3) + t.*Cov(1,4) + t.*Cov(2,3) + t.^2.*Cov(2,4);
%some necessary moments
% secondMu_m(1,:) = var_theta(1,:) + theta(1,:).^2;
% secondMu_d(1,:) = var_theta(2,:) + theta(2,:).^2;
secondMu_mMu_d(1,:) = covMu_mMu_d(1,:) + theta(1,:).*theta(2,:);
thirdMu_m2Mu_d(1,:) = theta(1,:).^2.*theta(2,:) + 2*covMu_mMu_d(1,:).*theta(1,:) + var_theta(1,:).*theta(2,:);
thirdMu_mMu_d2(1,:) = theta(2,:).^2.*theta(1,:) + 2*covMu_mMu_d(1,:).*theta(2,:) + var_theta(2,:).*theta(1,:);
fourthMu_m2Mu_d2(1,:) = theta(1,:).^2 .*theta(2,:).^2 + var_theta(1,:).*var_theta(2,:) + 2*covMu_mMu_d(1,:).^2 + ...
    4*covMu_mMu_d(1,:).*theta(1,:).*theta(2,:) + var_theta(1,:).*theta(2,:).^2 + var_theta(2,:).*theta(1,:).^2;

%the total methylation of hemimethylated sites (lambda)
theta(numOfRates+1,:) = theta(1,:) + theta(2,:) - theta(1,:).*theta(2,:);

%quantities refer to the appendix
varMu_mMu_d(1,:) = fourthMu_m2Mu_d2(1,:) - secondMu_mMu_d(1,:).^2;
CovMu_m2Mu_d(1,:) = thirdMu_m2Mu_d(1,:) - theta(1,:).*secondMu_mMu_d(1,:);
CovMu_mMu_d2(1,:) = thirdMu_mMu_d2(1,:) - theta(2,:).*secondMu_mMu_d(1,:);

%sum all the previous to get the variance of lambda
var_theta(numOfRates+1,:) = var_theta(1,:) + var_theta(2,:) + varMu_mMu_d(1,:) + 2*(covMu_mMu_d(1,:) - CovMu_m2Mu_d(1,:) - CovMu_mMu_d2(1,:));

%plot the standard deviation of lambda
errorbar(t, theta(numOfRates+1,:), sqrt(var_theta(numOfRates+1,:)), 'Marker', 'o', 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2(numOfRates+1,:));

% legend(labels, 'Location', 'east', 'FontSize', 12);
% legend('boxoff')
% set(legend5,...
%     'Position', [0.654281255142446 0.413078149920255 0.0453125 0.106858054226475],...   
%     'FontSize', 12);


%% Create axes4
axes4 = axes('Parent',figure1, ...
    'XTick', 1:numOfDays,...
    'XTickLabel', days, ...
    'FontSize', 12,...
    'Position', [0.183767361111111 0.102840909090909 0.08875 0.200715709728868]);
%'Position',[0.206736353077816 0.102840909090909 0.372822299651568 0.7]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes4,[0 1]);
box(axes4,'on');
hold(axes4,'on');


% Create ylabel
if strcmp(regionName(end), int2str(1))
    yLabelName = 'level per states';
else
    yLabelName = '';
end
ylabel(yLabelName,'FontSize', 15);

if (strcmp(dataType, 'oxBS'))
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
    bar1 = bar(plotedStates, 'BarLayout','stacked');
    set(bar1(1),'DisplayName','fullymethylated','FaceColor', [0.85 0.325 0.098]);
    set(bar1(2),'DisplayName','hemimethylated','FaceColor', [0.466 0.674 0.188]);
    set(bar1(3),'DisplayName','unmethylated','FaceColor', [0 0.447 0.741]);
    
end

% Create legend
% leg4 = legend(axes4, 'Location', 'best');
% set(leg4, 'FontSize', 10);


%% Create axes
axes5 = axes('Parent',figure1,'XTickLabel', dataLabels,...
    'XTick',1:numOfDays,...
    'XTickLabel', days, ...
    'FontSize', 12,...
    'Position',[0.29484375 0.102840909090909 0.06875 0.200715709728868]);


%Uncomment the following line to preserve the X-limits of the axes
xlim(axes5,[0 numOfDays+1]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes5,[0, 0.40]); %as in the paper
%ylim(axes5,[0, max(sum(hydroxyStates, 2)) + 0.05]);
box(axes5,'on');
hold(axes5,'on');

% Create ylabel
%ylabel('hydroxylation level', 'HorizontalAlignment', 'center', 'FontSize', 14);

if (strcmp(dataType, 'oxBS'))
    % Create multiple lines using matrix input to bar
    bar2 = bar(hydroxyStates,'Parent', axes5, 'BarLayout', 'stacked');
    set(bar2(1),'DisplayName','hm-mh');
    set(bar2(2),'DisplayName','uh-hu');
    set(bar2(3),'DisplayName','hh');

    % Create legend
%     legend1 = legend(axes5,'show');
%     set(legend1,...
%     'Position',[0.426462097167964 0.435833326351078 0.0688680013020834 0.0633333333333334],...
%     'FontSize',12);


end


end
