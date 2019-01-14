function varargout = UIhydroxy(varargin)
% UIHYDROXY MATLAB code for UIhydroxy.fig
%      UIHYDROXY, by itself, creates a new UIHYDROXY or raises the existing
%      singleton*.
%
%      H = UIHYDROXY returns the handle to a new UIHYDROXY or the handle to
%      the existing singleton*.
%
%      UIHYDROXY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UIHYDROXY.M with the given input arguments.
%
%      UIHYDROXY('Property','Value',...) creates a new UIHYDROXY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UIhydroxy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UIhydroxy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UIhydroxy

% Last Modified by GUIDE v2.5 28-Mar-2017 15:56:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UIhydroxy_OpeningFcn, ...
                   'gui_OutputFcn',  @UIhydroxy_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% --- Executes just before UIhydroxy is made visible.
function UIhydroxy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UIhydroxy (see VARARGIN)

dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'UIdialogBox'));
if (isempty(mainGuiInput)) ...
    || (length(varargin) <= mainGuiInput) ...
    || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.dialogBox = varargin{mainGuiInput+1};
    
    % Obtain handles of dialog window using GUIDATA with the caller's handle 
    dHandles = guidata(handles.dialogBox);
    
end

% % Update handles structure
% guidata(hObject, handles);

if dontOpen
   disp('-----------------------------------------------------');
   disp('Improper input arguments. Pass a property value pair') 
   disp('whose name is "changeme_main" and value is the handle')
   disp('to the changeme_main figure, e.g:');
   disp('   x = changeme_main()');
   disp('   changeme_dialog(''changeme_main'', x)');
   disp('-----------------------------------------------------');
% else
%    uiwait(hObject);
end


% UIWAIT makes UIhydroxy wait for user response (see UIRESUME)
% uiwait(handles.figure1);

figure(handles.figure1);
% set(handles.figure1, 'Visible', 'off');

%write the static text region
regionName = dHandles.FileName{1}(1:end-7);
set(handles.data_text1, 'String', regionName);
set(handles.data_text1, 'Visible', 'on');

%textbox in non oxBS case
% str = {'No oxBS data was given', 'for this region'};
% dim1 = [0.5 .5 .3 .3];  
% dim2 = [0.8 .08 .3 .3];
% handles.annot1 = annotation('textbox',dim1,'String',str,'FitBoxToText','on');
% handles.annot2 = annotation('textbox',dim2,'String',str,'FitBoxToText','on');
% set(handles.annot1, 'Visible', 'off');
% set(handles.annot2, 'Visible', 'off');

% handles.line_axes9 = axes;
% set(handles.line_axes9, 'Visible', 'off');
% set(handles.line_axes9, 'Position', [0.27, 0, 0.0000001, 1]);
% Switch off autoscaling.
% set(handles.line_axes9, 'Xlim', [0, 1], 'YLim', [0, 1]);
%create gui's line
% line([0, 0], [0, 1], 'Parent', handles.line_axes9);

%get handles.CpGFlags array
handles.CpGFlags = dHandles.CpGflags;

%get the CpGs with handles.CpGFlags == 1 and put them in the popup menu
CpGStrings{1} = 'Aggr. Data';
i=2;
for k=dHandles.CpGs
    if handles.CpGFlags(k) == 1
        CpGStrings{i} = strcat('CpG', num2str(k));
        i=i+1;
    end        
end        
handles.popup = uicontrol('Style', 'popup', 'String', CpGStrings, ...
    'Units', 'normalized', 'Position', [0.319 0.942 0.09 0.04], 'Callback', @popup_Callback);

%create the two panels for the single CpG plots and the scrollbar to work
panel1 = uipanel('Parent', handles.figure1);
panel2 = uipanel('Parent', panel1);

set(panel1, 'Units', 'normalized', 'Position', [0 0 0.27 1], 'Title', 'Single CpGs',...
    'FontSize', 14, 'BackgroundColor', [0.757, 0.867, 0.776]);

set(panel2, 'Units', 'normalized', 'Position', [0 -1 1 2], 'BackgroundColor', [0.757, 0.867, 0.776]);

%create scroll bar
s = uicontrol('Style', 'Slider', 'Parent', handles.figure1,...
      'Units','normalized','Position',[0.26 0 0.01 1], ...
      'BackgroundColor', [0.757, 0.867, 0.776], 'Value', 1, 'Min', 0, 'Max', 1, ... 
      'Callback', {@slider_callback1, panel2});


%% create the single CpGs plots

if any(handles.CpGFlags)

    txt = zeros(1, 14);
    for k=1:size(dHandles.CpGs, 2)
        cpg = dHandles.CpGs(k);
        if handles.CpGFlags(cpg) == 1
            % Add a text uicontrol to label CpG number.
            txt(cpg) = uicontrol('Style', 'text', 'parent', panel2, 'FontWeight', 'bold', 'Units', 'normalized',...
                'Position',[0 0.98-(cpg-1)*0.075 0.15 0.012], 'HorizontalAlignment', 'Left', 'BackgroundColor', [0.757, 0.867, 0.776],...
                'String', strcat('CpG', num2str(cpg)));
            
        days = dHandles.dataPoints;
        numOfDays = size(days, 1);

        pAllStates = dHandles.pAllStatesCpG(:,:,k);
        numOfExperiments = size(dHandles.pData, 3);
            
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

          
            
            for j=1:3
                if j==1
                    h = axes('Parent', panel2,  'Units', 'normalized', 'Position', [0.13, 0.927-(cpg-1)*0.075, 0.25 0.05]);
                    
                    %create xTickLabel
                    h.XTickLabel = days;
                    h.XTick = 1:numOfDays;

                                      
                    % Uncomment the following line to preserve the Y-limits of the axes
                    %Uncomment the following line to preserve the X-limits of the axes
                    xlim(h, [0 numOfDays+1]);
                    ylim(h, [0 1]);
                    box(h,'on');
                    hold(h,'on');

              
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
                        bar1 = bar(plotedStates, 'BarLayout', 'stacked');
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
                    
                    h.FontUnits = 'pixels';
                    h.FontSize = 8;

                elseif j==2
                    
                    h = axes('Parent', panel2, 'Units', 'normalized', 'Position', [0.44, 0.927-(cpg-1)*0.075, 0.2 0.05]);

                    h.XTickLabel = days;
                    h.XTick = 1:numOfDays;
                    

                    colormap(h, [0.929 0.5 0.125;0.929 0.694 0.125;0.929 0.9 0.125]);
                    %Uncomment the following line to preserve the X-limits of the axes
                    xlim(h, [0 numOfDays+1]);
                    % Uncomment the following line to preserve the Y-limits of the axes
                    %ylim(h,[0, max(sum(hydroxyStates, 2)) + 0.1]);
                    ylim(h, [0, 0.8]); 
                    box(h,'on');
                    hold(h,'on');

                    if numOfExperiments > 2
                        % Create multiple lines using matrix input to bar
                        bar2 = bar('Parent', h, [hydroxyStates, formalStates], 'BarLayout', 'stacked');
                        set(bar2(1),'DisplayName','hm-mh', 'Facecolor', [0.929 0.5 0.125]);
                        set(bar2(2),'DisplayName','uh-hu', 'Facecolor', [0.929 0.694 0.125]);
                        set(bar2(3),'DisplayName','hh', 'Facecolor', [0.929 0.9 0.125]);
                        set(bar2(4),'DisplayName','mf-fm', 'FaceColor', rgb('LightSkyBlue'));
                        set(bar2(5),'DisplayName','uf-fu', 'FaceColor', rgb('DeepSkyBlue'));
                        set(bar2(6),'DisplayName','hf-fh', 'FaceColor', rgb('DodgerBlue'));
                        set(bar2(7),'DisplayName','ff', 'FaceColor', rgb('CornflowerBlue'));
                    elseif numOfExperiments > 1
                        bar2 = bar('Parent', h, [hydroxyStates], 'BarLayout', 'stacked');
                        set(bar2(1),'DisplayName','hm-mh', 'Facecolor', [0.929 0.5 0.125]);
                        set(bar2(2),'DisplayName','uh-hu', 'Facecolor', [0.929 0.694 0.125]);
                        set(bar2(3),'DisplayName','hh', 'Facecolor', [0.929 0.9 0.125]);
                    else
                    %     set(dHandles.annot2, 'Visible', 'on');
                    end
                    
                    h.FontUnits = 'pixels';
                    h.FontSize = 8;
                 
                    
                else 
                    if  get(dHandles.popupmenu_eff, 'Value') ~= 2
                        
                        h = axes('Parent', panel2, 'Units', 'normalized', ...
                            'Position', [0.705, 0.927-(cpg-1)*0.075, 0.25 0.05]);


                        h.XTick = 1:days(end);

                        %names to variables
                        x = dHandles.xminmCpG(k,:);
                        Cov = dHandles.CovCpG(:,:,k);
                        t = days;
                        numOfDays = size(t,1);

                        %take out hydroxylation probability in case it is included
                        if strcmp(dHandles.paramNames{end}, 'p')
                            numOfParams = size(x, 2)-1;
                        else
                            numOfParams = size(x, 2);
                        end    

                        %find the degree (or numOfCoef) of the underlying polynom
                        numOfCoefPerRate = dHandles.maxDegree + 1;
                        % maxDegree = numOfCoef - 1;

                        numOfRates = floor(numOfParams/numOfCoefPerRate);
                        ratesNames = cell(1, numOfRates);
                        %get the rates names from paramNames vector
                        for i=1:numOfCoefPerRate:numOfParams
                            ratesNames{floor(i/numOfCoefPerRate) + 1} = dHandles.paramNames{i}(4:end);
                        end


                        rgbColor = [0.8500    0.3250    0.0980; %Red
                                 0         0.4470    0.7410;  %DarkBlue
                                 0.9290    0.6940    0.1250;  %yellow
                                 0         0.20    0.5410;      %lightblue
                                 rgb('Teal');                   %green
                                 0.5430    0         0   ];     %DarkRed


                        tPlot = 0:0.5:days(end);
                        for i=1:numOfRates
                            %get a vector with all the coefficients (b_0, b_1, ..) of rate i 
                            coef = x((i-1)*numOfCoefPerRate+1: (i-1)*numOfCoefPerRate+numOfCoefPerRate);
                            %compute the efficiency of rate i at every point of tPlot 
                            y = getEfficiency(coef, tPlot);

                            plot(h, tPlot, y, 'LineWidth', 2, 'LineStyle', '--', 'Color', rgbColor(i,:) );

    %                         var_b0 = Cov(i,i);
    %                         var_b1 = Cov(i+1,i+1);
    %                         cov = Cov(i,i+1);
    %                         
    %                         varOfy = var_b0 + var_b1*tPlot.^2 + 2*cov*tPlot;

    %                         [l,p] = boundedline(tPlot, y, 1.96*sqrt(varOfy), 'alpha', 'cmap', rgbColor(find(strcmp(ratesNames(i), allRatesArray)),:));
    %                         outlinebounds(l,p);

                            hold(h, 'on');
                        end


                        %the next is to be put as an argument above if we want specific Ticks
                        xlim(h, [t(1)-0.08, t(numOfDays)+0.08]);
                        ylim(h, [0, dHandles.ubR]);
                        box(h, 'on');
                        hold(h, 'on');


                        %define and plot each rate as a function of time
    %                     theta = zeros(numOfRates+1, numOfDays);
    %                     var_theta = zeros(numOfRates+1, numOfDays);
    % 
    %                     for i=1:2:2*numOfRates-1
    %                         theta((i+1)/2,:) = x(i) + x(i+1)*t;
    %                         %b0 should equal b1 for every parameter
    %                         %theta((i+1)/2,1) = x(i) + x(i+1)*1;
    %                         var_b0 = Cov(i,i);
    %                         var_b1 = Cov(i+1,i+1);
    %                         cov = Cov(i,i+1);
    %                         var_theta((i+1)/2,:) = var_b0 + var_b1*t.^2 + 2*cov*t;
    %                         errorbar(h, t, theta((i+1)/2,:), sqrt(var_theta((i+1)/2,:)), 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2((i+1)/2,:))
    %                         hold on;
    %                     end
    %                     %store and plot the lambda = mu1+mu2-mu1*mu2 over time
    %                     %theta(k,:) is the k-th entry of v =[mu_m, mu_d, eta] 
    %                     covMu_mMu_d(1,:) = Cov(1,3) + t.*Cov(1,4) + t.*Cov(2,3) + t.^2.*Cov(2,4);
    %                     %some necessary moments
    %                     % secondMu_m(1,:) = var_theta(1,:) + theta(1,:).^2;
    %                     % secondMu_d(1,:) = var_theta(2,:) + theta(2,:).^2;
    %                     secondMu_mMu_d(1,:) = covMu_mMu_d(1,:) + theta(1,:).*theta(2,:);
    %                     thirdMu_m2Mu_d(1,:) = theta(1,:).^2.*theta(2,:) + 2*covMu_mMu_d(1,:).*theta(1,:) + var_theta(1,:).*theta(2,:);
    %                     thirdMu_mMu_d2(1,:) = theta(2,:).^2.*theta(1,:) + 2*covMu_mMu_d(1,:).*theta(2,:) + var_theta(2,:).*theta(1,:);
    %                     fourthMu_m2Mu_d2(1,:) = theta(1,:).^2 .*theta(2,:).^2 + var_theta(1,:).*var_theta(2,:) + 2*covMu_mMu_d(1,:).^2 + ...
    %                         4*covMu_mMu_d(1,:).*theta(1,:).*theta(2,:) + var_theta(1,:).*theta(2,:).^2 + var_theta(2,:).*theta(1,:).^2;
    % 
    %                     %the total methylation of hemimethylated sites (lambda)
    %                     theta(numOfRates+1,:) = theta(1,:) + theta(2,:) - theta(1,:).*theta(2,:);
    % 
    %                     %quantities refer to the appendix
    %                     varMu_mMu_d(1,:) = fourthMu_m2Mu_d2(1,:) - secondMu_mMu_d(1,:).^2;
    %                     CovMu_m2Mu_d(1,:) = thirdMu_m2Mu_d(1,:) - theta(1,:).*secondMu_mMu_d(1,:);
    %                     CovMu_mMu_d2(1,:) = thirdMu_mMu_d2(1,:) - theta(2,:).*secondMu_mMu_d(1,:);
    % 
    %                     %sum all the previous to get the variance of lambda
    %                     var_theta(numOfRates+1,:) = var_theta(1,:) + var_theta(2,:) + varMu_mMu_d(1,:) + 2*(covMu_mMu_d(1,:) - CovMu_m2Mu_d(1,:) - CovMu_mMu_d2(1,:));
    % 
    %                     %plot the standard deviation of lambda
    %                     errorbar(h, t, theta(numOfRates+1,:), sqrt(var_theta(numOfRates+1,:)), 'LineStyle', '--', 'LineWidth', 2, 'Color', rgbcolor2(numOfRates+1,:));     

                        h.FontUnits = 'pixels';
                        h.FontSize = 8;
                    end 
                    
                end
                              
            end
        end
    end

end    


%% ----create the 5 main axes of gui----
%define an array of the 5 axes
axesArray = [handles.axes1, handles.axes2, handles.axes3, handles.axes4, handles.axes5, handles.axes6, handles.axes7];

%clear the axes
cla(handles.axes1, 'reset'); 
cla(handles.axes2, 'reset');
cla(handles.axes3, 'reset');
cla(handles.axes4, 'reset');
cla(handles.axes5, 'reset');
cla(handles.axes6, 'reset');
cla(handles.axes7, 'reset');
%clear legends
legend(handles.axes1,'hide');
legend(handles.axes2,'hide');
legend(handles.axes3,'hide');
legend(handles.axes4,'hide');
legend(handles.axes5,'hide');
legend(handles.axes6,'hide');
legend(handles.axes7,'hide');


%call plotOutput with this array
plotOutputToGui(axesArray, dHandles.xminm, dHandles.pData, dHandles.pModel, dHandles.pAllStates,...
    dHandles.CovPlot, dHandles.dataPointsName, dHandles.dataPoints, dHandles.labels, dHandles.paramNames, ...
    dHandles.knots, dHandles.maxDegree, dHandles.ubR);

%test call to confidenceEllipse
% plotConfRegions(dHandles.xminm, dHandles.CovPlot);


set(handles.figure1, 'Visible', 'on');

guidata(hObject, handles);




%% --- Executes on button press in exportpdf_pushbutton4.
function exportpdf_pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to exportpdf_pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the FigOut from the plotOutput
% figOut = plotOutput(handles.xminm, handles.pBis, handles.pOx, handles.pBisModel, handles.pOxModel, handles.pAllStates,...
%     handles.CovPlot, handles.dataPointsName, handles.dataPoints, handles.labels, handles.dataType);

% Obtain handles of dialog window using GUIDATA with the caller's handle 
dHandles = guidata(handles.dialogBox);

%get the name and the dir to store the .pdf file (default FileName is regionName.pdf)
%get the value of the popup window
str = get(handles.popup, 'String');
val = get(handles.popup,'Value');

%get the regionName for the outFileName
if strcmp(str{val}, 'Aggr. Data')
    %aggregated data
    paramsOut = dHandles.paramsOut;
    xminm = dHandles.xminm;
    CovPlot = dHandles.CovPlot;
    pData = dHandles.pData;
    pModel = dHandles.pModel;
    pAllStates = dHandles.pAllStates;
    
    regionName = dHandles.FileName{1}(1:end-7);
    
else    
    for k=dHandles.CpGs
        if handles.CpGFlags(k) == 1
            if strcmp(str{val}, strcat('CpG', num2str(k)))
                %CpGk data
                paramsOut = dHandles.paramsOutCpG(k,:);
                xminm = dHandles.xminmCpG(k,:);
                CovPlot = dHandles.CovPlotCpG(:,:,k);
                pData = dHandles.pDataCpG(:,:,k,:);
                pModel = dHandles.pModelCpG(:,:,k,:);
                %reshape pData, pModel for CPgk
                s = size(dHandles.pDataCpG(:,:,k,:));
                pData = reshape(pData, [s(1), s(2), s(4), s(3)]);                
                pModel = reshape(pModel, [s(1), s(2), s(4), s(3)]);
                
                pAllStates = dHandles.pAllStatesCpG(:,:,k);
                
                regionName = dHandles.FileName{1}(1:end-7); 
                regionName = strcat(regionName, '_CpG', num2str(k));
            end    
        end
    end
end

%get the dir to store the file
[handles.OutPdfFileName, handles.OutPdfPathName] = uiputfile('*.pdf', 'Save the .pdf output', strcat('~/Desktop/', regionName));

% (if no cancel was pressed)
if (handles.OutPdfFileName ~= 0) | (handles.OutPdfPathName  ~= 0)
    
    %picPath from uiputfile and regionName 
    picPath = strcat(handles.OutPdfPathName, handles.OutPdfFileName);
      
    fig = plotOutput(xminm, pData, pModel, pAllStates, CovPlot, dHandles.dataPointsName, dHandles.dataPoints, dHandles.labels, dHandles.paramNames, ...
        dHandles.knots, dHandles.maxDegree, dHandles.ubR);
    %export the figure as pdf in a3 format
    set(fig, 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape', 'PaperType', 'a3');
%     picPath = strcat('/Users/kyriakopou/Desktop/', regionName);
    print(fig, picPath, '-dpdf', '-r0')

    %call createPdfFile function to copy axes to pdf with path picPath 
%     createPdfFile(axesAndLegs, picPath, dHandles.dataType);
    
end


%% --- Executes on button press in exporttxt_pushbutton5.
function exporttxt_pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to exporttxt_pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain handles of dialog window using GUIDATA with the caller's handle 
dHandles = guidata(handles.dialogBox);

%get the value of the popup window
str = get(handles.popup, 'String');
val = get(handles.popup, 'Value');

%set the variables to the ones of the aggregated data
%or the corresponding CpG number
if strcmp(str{val}, 'Aggr. Data')
    paramsOut = dHandles.paramsOut;
    xminm = dHandles.xminm;
    Cov = dHandles.Cov;
    CovPlot = dHandles.CovPlot;
    sigma = dHandles.sigma;
    numOfDays = size(dHandles.dataPoints, 1);
    pData = dHandles.pData;
    pModel = dHandles.pModel;
    pAllStates = dHandles.pAllStates;
    apprCovFlag = dHandles.apprCovFlag;
    %regionName for the outputFileName
    regionName = dHandles.FileName{1}(1:end-7);
else
     %THIS CAN BE IMPROVED FOR SURE
    for k=1:size(dHandles.CpGs, 2)
        cpg = dHandles.CpGs(k);
        if handles.CpGFlags(cpg) == 1
            if strcmp(str{val}, strcat('CpG', num2str(cpg)))
                paramsOut = dHandles.paramsOutCpG(k,:);
                xminm = dHandles.xminmCpG(k,:);
                Cov = dHandles.CovCpG(:,:,k);
                CovPlot = dHandles.CovPlotCpG(:,:,k);
                sigma = dHandles.sigmaCpG(k,:);
                numOfDays = size(dHandles.dataPoints, 1);
                pData = dHandles.pDataCpG(:,:,k,:);
                pModel = dHandles.pModelCpG(:,:,k,:);
                %reshape pData, pModel for CPgk
                s = size(dHandles.pDataCpG(:,:,k,:));
                pData = reshape(pData, [s(1), s(2), s(4), s(3)]);                
                pModel = reshape(pModel, [s(1), s(2), s(4), s(3)]);
                
                pAllStates = dHandles.pAllStatesCpG(:,:,k);
                apprCovFlag = dHandles.apprCovFlagCpG(k);
                %regionName for the outputFileName
                regionName = dHandles.FileName{1}(1:end-7);
                regionName = strcat(regionName, '_CpG', num2str(cpg));
            end
        end    
    end
end

%get the name and the dir to store the .txt file (default FileName is regionName.txt)
% fileName = strcat(regionName, '_Results');
[handles.OutTxtFileName, handles.OutTxtPathName] = uiputfile('*.txt', 'Save the .txt output', strcat('~/Desktop/', regionName));

if (handles.OutTxtFileName ~= 0 | handles.OutTxtPathName ~= 0)
    
    outFilePath = strcat(handles.OutTxtPathName, handles.OutTxtFileName);
    fileID = fopen(outFilePath, 'wt', 'ieee-le','UTF-8');
    
    %create .txt file to fileID 
    createTxtFile(fileID, dHandles.paramNames, paramsOut, xminm, sigma, Cov, CovPlot, numOfDays, dHandles.labels, ... 
    pData, pModel, pAllStates, apprCovFlag);
   
    fclose(fileID);
end


guidata(hObject, handles);


%% -- popup callback function
function popup_Callback(hObject, eventdata, handles)

handles = guidata(hObject);

% Obtain handles of dialog window using GUIDATA with the caller's handle 
dHandles = guidata(handles.dialogBox);

str = get(hObject, 'String');
val = get(hObject,'Value');

axesArray = [handles.axes1, handles.axes2, handles.axes3, handles.axes4, handles.axes5, handles.axes6, handles.axes7];

%clear the axes
cla(handles.axes1, 'reset'); 
cla(handles.axes2, 'reset');
cla(handles.axes3, 'reset');
cla(handles.axes4, 'reset'); 
cla(handles.axes5, 'reset');
cla(handles.axes6, 'reset'); 
cla(handles.axes7, 'reset');
%clear legends
legend(handles.axes1,'hide');
legend(handles.axes2,'hide');
legend(handles.axes3,'hide');
legend(handles.axes4,'hide');
legend(handles.axes5,'hide');
legend(handles.axes6,'hide');
legend(handles.axes7,'hide');

%CONDITION IF IT CHANGED
% Set current data to the selected data set.
if strcmp(str{val}, 'Aggr. Data')
        
    plotOutputToGui(axesArray, dHandles.xminm, dHandles.pData, dHandles.pModel, dHandles.pAllStates,...
    dHandles.CovPlot, dHandles.dataPointsName, dHandles.dataPoints, dHandles.labels, dHandles.paramNames, dHandles.knots, dHandles.maxDegree, dHandles.ubR);
    
else
    %THIS CAN BE IMPROVED FOR SURE
    for k=1:size(dHandles.CpGs, 2)
        cpg = dHandles.CpGs(k);
        if handles.CpGFlags(cpg) == 1
            if strcmp(str{val}, strcat('CpG', num2str(cpg)))
                
                s = size(dHandles.pDataCpG(:,:,k,:));
                if length(s) > 3
                    pDataCpGk = reshape(dHandles.pDataCpG(:,:,k,:), [s(1), s(2), s(4), s(3)]);                
                    pModelCpGk = reshape(dHandles.pModelCpG(:,:,k,:), [s(1), s(2), s(4), s(3)]);
                else    
                    pDataCpGk = dHandles.pDataCpG(:,:,k,:);
                    pModelCpGk = dHandles.pModelCpG(:,:,k,:);
                end
                
                plotOutputToGui(axesArray, dHandles.xminmCpG(k,:), pDataCpGk, ...
                    pModelCpGk, dHandles.pAllStatesCpG(:,:,k), dHandles.CovPlotCpG(:,:,k), dHandles.dataPointsName, dHandles.dataPoints, ...
                    dHandles.labels, dHandles.paramNames, dHandles.knots, dHandles.maxDegree, dHandles.ubR);
            end
        end    
    end
end       
        
        
        
% Save the handles structure.
guidata(hObject,handles)


% --- Executes on button press in exportzip_pushbutton3.
function exportzip_pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to exportzip_pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to exporttxt_pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain handles of dialog window using GUIDATA with the caller's handle 
dHandles = guidata(handles.dialogBox);

%we have to run this for aggreg. data params values and for
%param values of each CpG


%get the directory to store the .zip file with the .txt and .pdf results
folder_name = uigetdir('~/Desktop', 'Select the folder to store the .zip file');

if (folder_name ~= 0)
    
    %aggregated data
    paramsOut = dHandles.paramsOut;
    xminm = dHandles.xminm;
    Cov = dHandles.Cov;
    CovPlot = dHandles.CovPlot;
    sigma = dHandles.sigma;
    paramNames = dHandles.paramNames;
    numOfDays = size(dHandles.dataPoints, 1);
    pData = dHandles.pData;
    pModel = dHandles.pModel;
    pAllStates = dHandles.pAllStates;
    apprCovFlag = dHandles.apprCovFlag;
    %regionName for the outputFileName
    regionName = dHandles.FileName{1}(1:end-7);
    
    folder_name = strcat(folder_name, '/', regionName, '_all');
    switch get(dHandles.popupmenu_eff, 'Value')
        case 1
            switch dHandles.model 
                case 1
                    folder_name = strcat(folder_name, '_eff_dtmc');
                case 2
                    folder_name = strcat(folder_name, '_eff_ctmc');
            end
        case 2
            folder_name = strcat(folder_name, '_levels');
    end    
    mkdir(folder_name);
    %get the name and the dir to store the .txt file (default FileName is regionName.txt)
    fileNameTxtAg = strcat(regionName, '.txt');
    
    %full path of the .txt file
    txtFilePath = strcat(folder_name, '/', fileNameTxtAg);
    fileID = fopen(txtFilePath, 'wt');
    
    createTxtFile(fileID, paramNames, paramsOut, xminm, sigma, Cov, CovPlot, numOfDays, dHandles.labels, ... 
            pData, pModel, pAllStates, apprCovFlag);
                 
    fclose(fileID);
    %get the figOut for the pdf
    figOut = plotOutput(xminm, pData, pModel, pAllStates, CovPlot, ...
    dHandles.dataPointsName, dHandles.dataPoints, dHandles.labels, dHandles.paramNames, dHandles.knots, dHandles.maxDegree, dHandles.ubR);
   
   
    %get fileNames of all files
    fileNamePdfAg = strcat(regionName, '.pdf');
    %full path of the .pdf file
    picPath = strcat(folder_name, '/', fileNamePdfAg);
    %open file with this path and prin the plots
    fileID = fopen(picPath, 'w');
    set(figOut, 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape', 'PaperType', 'a3', ...
        'PaperUnits', 'normalized');
    print(figOut, picPath, '-dpdf', '-r0');
    fclose(fileID);
            
    %do the same for all active CpGs    
    %THIS CAN BE IMPROVED FOR SURE
    fileNamesTxt = cell(1, size(dHandles.CpGs, 2));
    fileNamesPdf = cell(1, size(dHandles.CpGs, 2));
    
    for k=1:size(dHandles.CpGs, 2)
        cpg = dHandles.CpGs(k);
        if handles.CpGFlags(cpg) == 1
            paramsOut = dHandles.paramsOutCpG(k,:);
            xminm = dHandles.xminmCpG(k,:);
            Cov = dHandles.CovCpG(:,:,k);
            CovPlot = dHandles.CovPlotCpG(:,:,k);
            sigma = dHandles.sigmaCpG(k,:);
            paramNames = dHandles.paramNames;
            numOfDays = size(dHandles.dataPoints, 1);
%             pData = dHandles.pDataCpG(:,:,k,:);
%             pModel = dHandles.pModelCpG(:,:,k,:);
            pAllStates = dHandles.pAllStatesCpG(:,:,k);
            apprCovFlag = dHandles.apprCovFlagCpG(k);
            %regionName for the outputFileName
            regionName = dHandles.FileName{1}(1:end-7);
            regionNameCpG = strcat(regionName, '_CpG', num2str(cpg));
            %get the name and the dir to store the .txt file (default FileName is regionName.txt)
            fileNamesTxt{cpg} = strcat(regionNameCpG, '.txt');
        
            %full path of the .txt file
            txtFilePath = strcat(folder_name, '/', fileNamesTxt{cpg});
            fileID = fopen(txtFilePath, 'wt');
            
            s = size(dHandles.pDataCpG(:,:,k,:));
            pDataCpGk = reshape(dHandles.pDataCpG(:,:,k,:), [s(1), s(2), s(4), s(3)]);                
            pModelCpGk = reshape(dHandles.pModelCpG(:,:,k,:), [s(1), s(2), s(4), s(3)]);
            
            createTxtFile(fileID, paramNames, paramsOut, xminm, sigma, Cov, CovPlot, numOfDays, dHandles.labels, ... 
            pDataCpGk, pModelCpGk, pAllStates, apprCovFlag);
            
            fclose(fileID);
            
            %get the figOut for the pdf calling plotOutput
            figOut = plotOutput(xminm, pDataCpGk, pModelCpGk, pAllStates, CovPlot, ...
            dHandles.dataPointsName, dHandles.dataPoints, dHandles.labels, dHandles.paramNames, dHandles.knots, dHandles.maxDegree, dHandles.ubR);
 
            %get fileNames of all files
            fileNamesPdf{cpg} = strcat(regionNameCpG, '.pdf');
            %full path of the .pdf file
            picPath = strcat(folder_name, '/', fileNamesPdf{cpg});
            fileID = fopen(picPath, 'w');
            set(figOut, 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape', 'PaperType', 'a3',...
                'PaperUnits', 'normalized');
            print(figOut, picPath, '-dpdf', '-r0');
            fclose(fileID);        
        end    
    end
    
    %zip the created folder 
    zip(folder_name, folder_name);
    %remove the created folder
    rmdir(folder_name, 's');
    
    %get the cell array of all produced fileNames
%     allFileNames = [fileNamePdfAg, fileNameTxtAg, fileNamesPdf, fileNamesTxt];
%     allFileNames = allFileNames(~cellfun('isempty', allFileNames));   
    %create the zip file and put all results in there    
%     zipFileFullName = strcat(folder_name, '/', dHandles.FileName{1}(1:end-7));
    %delete the produced files
%     allFilesPath = strcat(folder_name, '/', allFileNames);
%     for i=1:size(allFilesPath, 2)
%         delete(allFilesPath{i});
%     end
    %delete initial files
%     txtFilesPath = strcat(folder_name, '/', '*.txt');
%     pdfFilesPath = strcat(folder_name, '/', '*pdf');

    
end


guidata(hObject, handles);



function slider_callback1(src, eventdata, arg1)
val = get(src,'Value');
set(arg1,'Position', [0 -val 1 2]);  


function createTxtFile(fileID, paramNames, paramsOut, xminm, sigma, Cov, CovPlot, numOfDays, labels, ... 
    pData, pModel, pAllStates, apprCovFlag)
    
    numOfExperiments = sum(any(any(pData))); 
        
    if ~isempty(xminm)
        
        %Print the estimated parameters the stds
        %and relative half widths
        fprintf(fileID, 'The parameters of the model are: \n');
        fprintf(fileID, 'l: Total methylation on hemi-methylated CpGs \n');
        fprintf(fileID, 'm: Maintenance methylation \n');
        fprintf(fileID, 'd: De novo methylation \n');
        if numOfExperiments > 1 
            fprintf(fileID, 'e: Hydroxylation \n');
            fprintf(fileID, 'p: Prob. that Dnmt1 does not recognize 5hmC \n\n');
        end
        dataFormatSpec = ['%s = %4.3f ' char(177) ' (%4.3f)\n'];
        for i = 1:size(paramsOut, 2)
            %fprintf(fileID, formatSpec, paramNames{1,:});
            %formatSpec = '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f\n';
            fprintf(fileID, dataFormatSpec, char(paramNames(i)), paramsOut(i), sigma(i));
        end    
        fprintf(fileID, '\n');

        % fprintf('The relative half widths are: \n');
        % relHalfWidth = abs(1.96*sigma ./ params);
        % display(relHalfWidth);
        % fprintf('\n');

        %compute the lambda coefficients (need all params and all stds)
        lambdaCoef(xminm, CovPlot, fileID)


        if (numOfDays > 2)
            if numOfExperiments > 1
                numOfStates = 16;
                %define the g function and the Jacobian for the rate params of interest (not-nuisance)
                Wald = zeros(1, floor(size(paramsOut, 2) / 2));
                for testParam = 2:2:size(paramsOut, 2)
                    g = paramsOut(testParam);
                    Jacob = zeros(1, size(paramsOut,2));
                    Jacob(testParam) = 1;
    %                 Cov = nan(size(testParam)); %only for the stuff now to work
                    Wald(testParam / 2) = g * inv(Jacob * Cov * Jacob')*g';
                    %convert cell array position to string (char array)
                    paramName = char(paramNames(testParam));
                    fprintf(fileID, 'The Wald statistic for the %s parameter is: %d \n', paramName, Wald(testParam / 2));
                end

                %wald test for non-recognition prob. if it is in the keeped params
    %             if (strcmp(paramNames{end}, 'p') )
    %                 gP = paramsOut(end) - 1;
    %                 JacobP = zeros(1, size(paramsOut, 2));
    %                 JacobP(7) = 1;
    %                 WaldP = gP * inv(JacobP*Cov*JacobP')*gP';
    %                 fprintf(fileID, 'The Wald statistic for the p parameter is: %d\n', WaldP);
    %             end

                %wald test for total methylation lambda
                gLambda = [paramsOut(2) + paramsOut(4) - paramsOut(2)*paramsOut(3) - paramsOut(1)*paramsOut(4);
                            paramsOut(2)*paramsOut(4)];
                JacobLambda = [-paramsOut(4), 1-paramsOut(3), -paramsOut(2), 1-paramsOut(1), zeros(1, size(paramsOut, 2)-4);
                               0, paramsOut(4), 0, paramsOut(2), zeros(1, size(paramsOut, 2) - 4)];
                WaldLambda = gLambda' * inv(JacobLambda*Cov*JacobLambda')*gLambda;
                fprintf(fileID, 'The Wald statistic for total methylation is: %d\n', WaldLambda);

            else

                %define the g function and the Jacobian for the rate params of
                %interest (here it is only b1_maintenance and b1_denovo )
                for testParam = 2:2:4
                    g = paramsOut(testParam);
                    Jacob = zeros(1, size(paramsOut,2));
                    Jacob(testParam) = 1;
                    Wald = g * inv(Jacob * Cov * Jacob')*g';
                    %convert cell array position to string (char array)
                    paramName = char(paramNames(testParam));
                    fprintf(fileID, 'The Wald statistic for the %s parameter is: %d \n', paramName, Wald);
                end

                %wald test for total methylation lambda
                gLambda = [paramsOut(2) + paramsOut(4) - paramsOut(2)*paramsOut(3) - paramsOut(1)*paramsOut(4);
                            paramsOut(2)*paramsOut(4)];
                JacobLambda = [-paramsOut(4), 1-paramsOut(3), -paramsOut(2), 1-paramsOut(1), zeros(1, size(paramsOut, 2)-4);
                               0, paramsOut(4), 0, paramsOut(2), zeros(1, size(paramsOut, 2) - 4)];
                WaldLambda = gLambda' * inv(JacobLambda*Cov*JacobLambda')*gLambda;
                fprintf(fileID, 'The Wald statistic for total methylation is: %d\n', WaldLambda);

            end
        end

        if apprCovFlag
            fprintf(fileID, '\n');
            fprintf(fileID, 'Warning! stds are approximations after computing \nthe nearest symmetric-positive definite Cov matrix \n');                
        end

        fprintf(fileID, '---------------------------------------------------------------\n');
    
    end
    

    C = num2cell(pData(:,:,1));
    C = [labels, C];
    dataFormatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';
    
    fprintf(fileID, 'The data distribution of BS states \n\n');
    fprintf(fileID, 'time TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, dataFormatSpec,C{row,:});
    end
    fprintf(fileID, '\n');
    
    C = num2cell(pModel(:,:,1));
    C = [labels, C];
    dataFormatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';
    
    fprintf(fileID, 'The predicted distribution of BS states \n\n');
    fprintf(fileID, 'time TT TC CT CC \n');
    [nrows, ncols] = size(C);
    for row = 1:nrows
        fprintf(fileID, dataFormatSpec,C{row,:});
    end
    
    fprintf(fileID, '---------------------------------------------------------------\n');
    if numOfExperiments > 1
        numOfStates = 16;
        C = num2cell(pData(:,:,2));
        C = [labels, C];
        dataFormatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The data distribution of oxBS states \n\n');
        fprintf(fileID, 'time TT TC CT CC \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, dataFormatSpec,C{row,:});
        end
        fprintf(fileID, '\n');

        C = num2cell(pModel(:,:,2));
        C = [labels, C];
        dataFormatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The predicted distribution of oxBS states \n\n');
        fprintf(fileID, 'time TT TC CT CC \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, dataFormatSpec,C{row,:});
        end
    
        fprintf(fileID, '---------------------------------------------------------------\n');
    end
    
    
    if numOfExperiments > 2
        numOfStates = 16;
        C = num2cell(pData(:,:,3));
        C = [labels, C];
        dataFormatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The data distribution of mabBS states \n\n');
        fprintf(fileID, 'time TT TC CT CC \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, dataFormatSpec, C{row,:});
        end
        fprintf(fileID, '\n');

        C = num2cell(pModel(:,:,2));
        C = [labels, C];
        dataFormatSpec = '%s %1.4f %1.4f %1.4f %1.4f\n';

        fprintf(fileID, 'The predicted distribution of mabBS states \n\n');
        fprintf(fileID, 'time TT TC CT CC \n');
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fileID, dataFormatSpec,C{row,:});
        end
    
        fprintf(fileID, '---------------------------------------------------------------\n');
    end
    
    C = num2cell(pAllStates);
    C = [labels, C];
    
    headerFormatSpec = [repmat('%s ', [1, numOfStates+1]), '\n'];
    dataFormatSpec = ['%s ', repmat('%1.4f ', [1, numOfStates]), '\n'];
    
    stateSet = createStateSet({'u', 'm', 'h', 'f'}, 2);
    
    fprintf(fileID, 'The predicted distribution of the hidden states \n\n');
    header = ['time', stateSet];
    fprintf(fileID, headerFormatSpec, header{:});
    [nrows, ncols] = size(C);
    
    for row = 1:nrows
        fprintf(fileID, dataFormatSpec, C{row,:});
    end
    fprintf(fileID, '\n');
    
    
function createPdfFile(axesAndLegs, picPath, dataType)


% Create new figure FigOut with pdf appropriate size and copy everything there!
figOut = figure('Visible', 'off', 'Name', 'Model prediction', 'Units', 'Points', ...
    'Colormap',[0.929 0.5 0.125;0.929 0.694 0.125;0.929 0.9 0.125], ...
    'Position', [0, 0, 1152, 600]);

% figOut.Units = 'points';
% figOutPositionPts = figOut.Position
    
h(1) = copyobj(axesAndLegs(1), figOut);  
%put the column legend in the first axes
labels = {'TT', 'TT pr', 'TC', 'TC pr', 'CT', 'CT pr', 'CC', 'CC pr'};
% set(figOut, 'CurrentAxes', h(1));
% columnlegend(2, labels, 'fontsize', 7);
legendflex(labels, 'anchor', {'nw','nw'}, ...
        'buffer', [-5 -5], ...
        'ncol', 4, ...
        'fontsize', 6, ...
        'xscale', 0.6, ...
        'box', 'off');

%THIS IS HARDCODED
h(2) = copyobj(axesAndLegs(2), figOut);
h(3:4) = copyobj(axesAndLegs(3:4), figOut);
h(5:6) = copyobj(axesAndLegs(5:6), figOut);
if (strcmp(dataType, 'oxBS'))
    h(7:8) = copyobj(axesAndLegs(7:8), figOut);
end    
    
%set the position in the A4 pdf paper
set(h(1), 'Units', 'normalized', 'Position', [0.183767361111111 0.59006379585327...
    0.177951388888889 0.281180223285486]);  
set(h(2), 'Units', 'normalized', 'Position', [0.3900 0.59006379585327...
    0.17951388888889 0.281180223285486]);
set(h(4), 'Units', 'normalized', 'Position', [0.62084375 0.59006379585327...
    0.229582849399922 0.281180223285485]);
set(h(6), 'Units', 'normalized', 'Position', [0.183767361111111 0.102840909090909...
    0.38675 0.410715709728868]);
if (strcmp(dataType, 'oxBS'))
    set(h(8), 'Units', 'normalized', 'Position', [0.62084375 0.0994318181818182...
        0.229582849399922 0.410715709728868]);
    %increase ylabel font size for hydroxylation plot
    fontSizePts = ceil(get(h(8), 'FontSize'));
    ylabel(h(8), 'hydroxylation level', 'HorizontalAlignment', 'center', 'FontSize', ceil(fontSizePts*1.5));
end 


set(figOut, 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape', 'PaperType', 'a4', 'PaperUnits', 'normalized');
print(figOut, picPath, '-dpdf', '-r0');



%% --- Outputs from this function are returned to the command line.
function varargout = UIhydroxy_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = [];

% 
 %% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%


% --- Executes when user attempts to close figure.
function figure1_CloseRequestFcn(hObject,eventdata,handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UIhydroxy('UIdialogBox', handles.figure1);

% The GUI is no longer waiting, so destroy it now.
delete(hObject);
