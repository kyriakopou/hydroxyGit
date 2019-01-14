function varargout = UIdialogBox(varargin)
% UIDIALOGBOX MATLAB code for UIdialogBox.fig
%      UIDIALOGBOX, by itself, creates a new UIDIALOGBOX or raises the existing
%      singleton*.
%
%      H = UIDIALOGBOX returns the handle to a new UIDIALOGBOX or the handle to
%      the existing singleton*.
%
%      UIDIALOGBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UIDIALOGBOX.M with the given input arguments.
%
%      UIDIALOGBOX('Property','Value',...) creates a new UIDIALOGBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UIdialogBox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UIdialogBox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UIdialogBox

% Last Modified by GUIDE v2.5 29-May-2017 11:07:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UIdialogBox_OpeningFcn, ...
                   'gui_OutputFcn',  @UIdialogBox_OutputFcn, ...
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


%% --- Executes just before UIdialogBox is made visible.
function UIdialogBox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UIdialogBox (see VARARGIN)

% Choose default command line output for UIdialogBox
handles.output = hObject;

set(handles.params_pushbutton3, 'Enable','off');
% set(handles.exportpdf_pushbutton4, 'Enable','off');
% set(handles.exporttxt_pushbutton5, 'Enable','off');

% set(handles.data_text1, 'Visible', 'off');
set(handles.BS_text1, 'Visible', 'off');
set(handles.oxBS_text2, 'Visible', 'off');
set(handles.mabBS_text, 'Visible', 'off');
set(handles.errors_text3, 'Visible', 'off');

set(handles.tick_axes1, 'Visible', 'off');
set(handles.tick_axes2, 'Visible', 'off');
set(handles.tick_axes3, 'Visible', 'off');
set(handles.tick_axes4, 'Visible', 'off');


%read and convert example image 
imgQ = imread('questionmark.png', 'BackgroundColor',[1 1 1]);
imgQ = imresize(imgQ, 1/16);

%initialize undo image
[img, map, alpha] = imread('undo.png', 'BackgroundColor', [1 1 1]);
img = imresize(img, 1/32);
img = double(img)/255;
index1 = img(:,:,1) == 1;
index2 = img(:,:,2) == 1;
index3 = img(:,:,3) == 1;
indexWhite = index1+index2+index3==3;
for idx = 1 : 3
   rgb = img(:,:,idx);     % extract part of the image
   rgb(indexWhite) = NaN;  % set the white portion of the image to NaN
   img(:,:,idx) = rgb;     % substitute the update values
end

%initialize example buttons
set(handles.exampleBS_pushbutton7, 'CData', imgQ);
set(handles.exampleOxBS_pushbutton8, 'CData', imgQ);
set(handles.exampleMabBS_pushbutton9, 'CData', imgQ);
set(handles.errorsEx_pushbutton10, 'CData', imgQ);
%initialize undo button
set(handles.undo_pushbutton6, 'Enable', 'off', 'Visible', 'off', 'CData', img);
set(handles.undo_pushbutton7, 'Enable', 'off', 'Visible', 'off', 'CData', img);
%initialize FileName, dataType and 1st legend
handles.FileName = cell(1,3);
handles.errorsFileName = '';
handles.dataType = '';
%initialize CpGflags
handles.CpGflags = zeros(1,14);
%initialize old CpG number
handles.oldCpGs = [];
%regionNames
handles.regionNameBS = [];
handles.regionNameOx = [];
handles.regionNameMab = [];
%initialize the following fields to empty
handles.xminm = [];
handles.paramsOut = [];
handles.sigma = [];
handles.Cov = [];
handles.pData = [];
handles.pModel = [];
handles.pAllStates = [];
handles.Cov = [];
handles.CovPlot = [];
handles.apprCovFlag = [];

%TO BE CHECKED
% handles.paramNames = {'b0_m', 'b1_m', 'b0_d', 'b1_d', 'b0_e', 'b1_e', 'b0_f', 'b1_f', 'b0_dm', 'b1_dm', 'p'};
handles.modelFileName = 'modelActiveLin.txt';   
% handles.modelFileName = 'modelActiveQuad.txt';

%create Panel with dynamic size and dynamic 
%number of CpGs.
% handles.panel = uipanel('parent', UIdialogBox,...
% 'Title','single CpGs',... 
% 'position', [.113 .346 .80 .364]);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UIdialogBox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% --- Outputs from this function are returned to the command line.
function varargout = UIdialogBox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on button press in load_pushbutton1.
function load_pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.FileName{1}, handles.PathName] = uigetfile('*.txt', 'Select BS data', ...
    '/Users/kyriakopou/Documents/code/MATLAB/HydroxyMethylation/dataFilesNew/newFormat/' ,'Multiselect', 'off');

if (handles.FileName{1} ~= 0) %if no cancel was pressed
    
    %check the format of fileName (_BS.txt)
    if (~strcmp(handles.FileName{1}(end-6:end), '_BS.txt'))
        msg = 'Not appropriate input fileName';
        msgbox('Please give appropriate _BS.txt fileName');
        error(msg); 
    end
    
    %set the names of BS file
    file1 = strcat(handles.PathName, handles.FileName{1});
    
    %get the regionName from the fileName (difference btw two is that regionName is never 0)
    handles.regionNameBS = handles.FileName{1}(1:end-4);
    
    
    %if _oxBS file already given check if #rows of BS agrees with #rows of oxBS
    %and if not halt and return error msg -> number of data Points differs
    if (~isempty(handles.regionNameOx))
%         if (size(handles.dataPointsOx, 1) ~= size(handles.dataPoints, 1)) 
%             msg = 'Not appropriate input files';
%             msgbox('Number of dataPoints differs in BS and oxBS files');
%             error(msg); 
%         end
        %if _BS and _oxBS do not refer to the same region
        if (~strcmp(handles.regionNameBS(1:end-3), handles.regionNameOx(1:end-5)))
            msg = 'Not appropriate input fileName';
            msgbox({'The _BS.txt and ox_BS.txt files'; 'do not refer to the same locus.'});
            error(msg); 
        end         
    end
    
    %call readData 
%     [handles.dataPointsName, handles.dataPoints, handles.labels, handles.obsBis, ...
%         handles.errorsBis, handles.process] = UIreadData(file1);
    
    handles.obsData = [];
    handles.obsDataCpG = [];
    handles.CpGs = [];
    %call readDataNew
    [handles.dataPointsName, handles.dataPoints, ...
        handles.obsData(:,:,1), handles.obsDataCpG(:,:,:,1), handles.CpGs] = UIreadDataNew(file1);
    
    %initialize obsOx, obsOxCpG and errorsBis errorsOx in case we haven't
    %loaded oxBS and MAB_BS - we need this if bec we cant initiliaze before
    %since we don't know the size of obsMatrix
    if (isempty(handles.regionNameOx))
        handles.FileName{2} = '';
        handles.obsData(:,:,2) = zeros(size(handles.obsData(:,:,1)));
        handles.obsDataCpG(:,:,:,2) = zeros(size(handles.obsDataCpG(:,:,:,1)));
        handles.obsData(:,:,3) = zeros(size(handles.obsData(:,:,1)));
    end
    
    %get the number of days
    numOfDays = size(handles.dataPoints, 1);
    
    %if no errors file has been given labels = dataPointsName + 0 | 1 |...
    %process = rep, ..., rep and errors are set to 0.005 and 0.07
    if isempty(handles.errorsFileName)
        %initialize errors 
        %Bis
        handles.errors = zeros(size(handles.dataPoints, 1), 3, 3);
        handles.errors(:,1,1) = 0.005;
        handles.errors(:,2:3,1) = 0.07;
        handles.errors(:,4,1) = 0.03;
        %oxBis
        handles.errors(:,1,2) = 0.005;
        handles.errors(:,2:3,2) = 0.07;
        handles.errors(:,4,3) = 0.03;
        %mabBis
        handles.errors(:,1,3) = 0.05;
        handles.errors(:,2:3,3) = 0.07;
        handles.errors(:,4,3) = 0.03;
        %initialize labels taking only the first letter (day ->d hour ->h)
        handles.labels = cellstr(repmat(handles.dataPointsName, numOfDays, 1));
        for i=1:numOfDays
            handles.labels{i} = strcat(handles.labels{i}(1), num2str(handles.dataPoints(i))); 
        end
        %initialize process
        a = {'rep'};
        handles.process = repmat(a, numOfDays, 1);
    end
    
    %switch on BS_text
    set(handles.BS_text1, 'String', handles.regionNameBS);
    set(handles.BS_text1, 'Visible', 'on');
    

    %put the tick image
    axes(handles.tick_axes2);
    Img = imread('tick.png', 'BackgroundColor', [0 0 0]);
    mask = bsxfun(@eq,Img,reshape([0 0 0], 1,1,3));
    handles.tick1 = image(Img,'alphadata', 1-double(all(mask,3)));
    
    axis off;
    axis image;  
        
    
    %put select all check box on the panel
    handles.selAll = uicontrol('parent', handles.uipanel, 'style', 'checkbox', ...
                   'tag', 'selAll', ...
                   'string', 'Select All',...
                   'Units','normalized',...
                   'position', [0.45, 0.03, 0.597, 0.108], ...
                   'callback', @selectAll_Callback);
    
    %delete and clear all old CpG checkboxes           
    for i = handles.oldCpGs;
        delete(handles.checkbox(i));
        handles.CpGflags(i) = 0;
    end           
    %put one check box for each CpG (new) of the data
    for i = handles.CpGs
        if i <= 7
            x = 0.05;
        else    
            x = 0.5;
        end
        
        position = [x, 0.8-0.105*mod(i-1,7), 0.40, 0.108];
        handles.checkbox(i) = uicontrol('parent', handles.uipanel, 'style', 'checkbox', ...
                   'tag', sprintf('checkbox%d', i), ...
                   'string', sprintf('CpG%d', i), ...
                   'Units','normalized',...
                   'position', position,...
                   'callback', str2func(sprintf('checkbox%d_Callback', i)));
    end               
    
    %update the CpG vector
    handles.oldCpGs = handles.CpGs;

    
    %switch off non-oxBS annotation text for oxBS data
%     set(handles.annot1, 'Visible', 'off');
%     set(handles.annot2, 'Visible', 'off');
    %make old text visible
%     set(handles.data_text1, 'Visible', 'off');
    %enable optimize button
    set(handles.params_pushbutton3, 'Enable','on'); 
    
    %disable output buttons
%     set(handles.exportpdf_pushbutton4, 'Enable', 'off');
%     set(handles.exporttxt_pushbutton5, 'Enable', 'off');
    
    
    
%     %print message box
%     msg = strcat(handles.FileName{1}, ' has been succesfully loaded!');
%     msgbox(msg);
    
end



guidata(hObject, handles);



%% --- Executes on button press in load_pushbutton2.
function load_pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[handles.FileName{2}, handles.PathName] = uigetfile('*.txt', 'Select oxBS data', ...
    '/Users/kyriakopou/Documents/code/MATLAB/HydroxyMethylation/dataFilesNew/newFormat/' ,'Multiselect', 'off');

%if no cancel was pressed
if (handles.FileName{2} ~= 0)

    %check that fileName ends with (_oxBS.txt)
    if (~strcmp(handles.FileName{2}(end-8:end), '_oxBS.txt'))
        msg = 'Not appropriate input fileName';
        msgbox('Please give appropriate _oxBS.txt fileName.');
        error(msg); 
    end
    
    
    %set the name of oxBS file
    file2 = strcat(handles.PathName, handles.FileName{2});

%     %call readData 
%     [handles.dataPointsName, handles.dataPointsOx, handles.labels,...
%         handles.obsOx, handles.errorsOx, handles.process] = UIreadData(file2);

    %get regionName
    handles.regionNameOx = handles.FileName{2}(1:end-4);
    
    
    %if BS file already given check if #rows of BS agrees with #rows of oxBS
    %and if not halt and return error msg -> number of data Points differs
    %if (~isempty(handles.FileName{1}) & handles.FileName{1} ~= 0)
    if (~isempty(handles.regionNameBS))    
%         if (size(handles.dataPointsOx, 1) ~= size(handles.dataPoints)) 
%             msg = 'Not appropriate input files';
%             msgbox('Number of dataPoints differs in BS and oxBS files');
%             error(msg); 
%         end
        
        %check that region is the same with _BS if _BS exists
        if (~strcmp(handles.regionNameBS(1:end-3), handles.regionNameOx(1:end-5)))
            msg = 'Not appropriate input fileName';
            msgbox({'The _BS.txt and ox_BS.txt files'; 'do not refer to the same locus.'});
            error(msg); 
        end
    end
    
    %call readDataNew 
    [handles.dataPointsName, handles.dataPointsOx, handles.obsData(:,:,2), handles.obsDataCpG(:,:,:,2)] = UIreadDataNew(file2);
    
    %write the oxBS txt box and make it visible
    %switch on BS_text
    set(handles.oxBS_text2, 'String', handles.regionNameOx);
    set(handles.oxBS_text2, 'Visible', 'on');
    
    %change the dataType to oxBS
    handles.dataType = 'oxBS';
    
    %if BS data already given enable params button
    if (~isempty(handles.regionNameBS))
        set(handles.params_pushbutton3, 'Enable','on');
    end
    
%     set(handles.annot1, 'Visible', 'off');
%     set(handles.annot2, 'Visible', 'off');
    
    %clear the axes
%     cla(handles.tick_axes2, 'reset');
%     cla(handles.tick_axes3, 'reset'); 
%     cla(handles.axes3, 'reset'); 
%     cla(handles.tick_axes1, 'reset'); 
%     cla(handles.axes5, 'reset');

%     %clear legends
%     legend(handles.tick_axes2,'hide');
%     legend(handles.tick_axes3,'hide');
%     legend(handles.axes3,'hide');
%     legend(handles.tick_axes1,'hide');
%     legend(handles.axes5,'hide');
    
    %enable undo button
    set(handles.undo_pushbutton6, 'Enable', 'on', 'Visible', 'on');
    %put the tick image
    axes(handles.tick_axes3);
    Img = imread('tick.png', 'BackgroundColor', [0 0 0]);
    mask = bsxfun(@eq,Img,reshape([0 0 0], 1,1,3));
    handles.tick2 = image(Img,'alphadata',1-double(all(mask,3)));
    %image(Img);
    axis off;
    axis image;  
   
    
%     %print message box
%     msg = strcat(handles.FileName{2}, ' has been succesfully loaded!');
%     msgbox(msg);
    
end
    
guidata(hObject, handles);


%% --- Executes on button press in load_pushbutton3.
function load_pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.FileName{3}, handles.PathName] = uigetfile('*.txt', 'Select mabBS data', ...
    '/Users/kyriakopou/Documents/code/MATLAB/HydroxyMethylation/dataFilesNew/newFormat/' ,'Multiselect', 'off');

%if no cancel was pressed
if (handles.FileName{3} ~= 0)

    %check that fileName ends with (_oxBS.txt)
    if (~strcmp(handles.FileName{3}(end-9:end), '_mabBS.txt'))
        msg = 'Not appropriate input fileName';
        msgbox('Please give appropriate _mabBS.txt fileName.');
        error(msg); 
    end
    
    
    %set the name of mabBS file
    file3 = strcat(handles.PathName, handles.FileName{3});

%     %call readData 
%     [handles.dataPointsName, handles.dataPointsOx, handles.labels,...
%         handles.obsOx, handles.errorsOx, handles.process] = UIreadData(file2);

    %get regionName
    handles.regionNameMab = handles.FileName{3}(1:end-4);
    
    
    %if BS or oxBS file already given check if #rows of BS or oxBS agrees with #rows of mabBS
    %and if not halt and return error msg -> number of data Points differs
    %if (~isempty(handles.FileName{1}) & handles.FileName{1} ~= 0)
    if (~isempty(handles.regionNameBS))    
%         if (size(handles.dataPointsOx, 1) ~= size(handles.dataPoints)) 
%             msg = 'Not appropriate input files';
%             msgbox('Number of dataPoints differs in BS and oxBS files');
%             error(msg); 
%         end
        
        %check that region is the same with _BS if _BS exists
        if (~strcmp(handles.regionNameBS(1:end-3), handles.regionNameMab(1:end-6)))
            msg = 'Not appropriate input fileName';
            msgbox({'The _BS.txt and ox_BS.txt files'; 'do not refer to the same locus.'});
            error(msg); 
        end
    end
    
    if (~isempty(handles.regionNameOx))    
%         if (size(handles.dataPointsOx, 1) ~= size(handles.dataPoints)) 
%             msg = 'Not appropriate input files';
%             msgbox('Number of dataPoints differs in BS and oxBS files');
%             error(msg); 
%         end
        
        %check that region is the same with _BS if _BS exists
        if (~strcmp(handles.regionNameOx(1:end-5), handles.regionNameMab(1:end-6)))
            msg = 'Not appropriate input fileName';
            msgbox({'The ox_BS.txt and mab_BS.txt files'; 'do not refer to the same locus.'});
            error(msg); 
        end
    end
    
    %call readDataNew 
    [handles.dataPointsName, handles.dataPointsMab, handles.obsData(:,:,3), handles.obsDataCpG(:,:,:,3)] = UIreadDataNew(file3);
    
    %write the oxBS txt box and make it visible
    %switch on BS_text
    set(handles.mabBS_text, 'String', handles.regionNameMab);
    set(handles.mabBS_text, 'Visible', 'on');
    
    %change the dataType to oxBS
    handles.dataType = 'mabBS';
    
    %if BS and oxBS data already given enable params button
    if (~isempty(handles.regionNameBS) && ~isempty(handles.regionNameOx))
        set(handles.params_pushbutton3, 'Enable','on');
    end
    
    
    %enable undo button
    set(handles.undo_pushbutton7, 'Enable', 'on', 'Visible', 'on');
    %put the tick image
    axes(handles.tick_axes4);
    Img = imread('tick.png', 'BackgroundColor', [0 0 0]);
    mask = bsxfun(@eq,Img,reshape([0 0 0], 1,1,3));
    handles.tick2 = image(Img,'alphadata',1-double(all(mask,3)));
    %image(Img);
    axis off;
    axis image;  
   
    
%     %print message box
%     msg = strcat(handles.FileName{2}, ' has been succesfully loaded!');
%     msgbox(msg);
    
end
    
guidata(hObject, handles);


%% --- Executes on button press in params_pushbutton3.
function params_pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to params_pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hf=findobj('Name', 'H(O)TA Main Window');
close(hf);

%get the selectedObject from the butting group
% blah = get(handles.opt_uibuttongroup1, 'SelectedObject');
% dataType = get(blah, 'String');

%get the dataType depending on whether oxBS was loaded or not
% if (~isempty(handles.FileName{2}))
%     handles.dataType = 'oxBS';
% end


%disable paramsbutton, enable plot and export buttons
set(handles.params_pushbutton3, 'Enable','off');
% set(handles.exportpdf_pushbutton4, 'Enable','on');
% set(handles.exporttxt_pushbutton5, 'Enable','on');

%message box for optimization
m = msgbox('Wait until optimization is done ...');
edithandle = findobj(m,'Style','pushbutton');
set(edithandle,'Visible','off');

tic; 

%call UIestimateDSHydroxy to get the optimal params and the model
%prediction for the aggregated data and the single CpGs checked
regionName = handles.FileName{1}(1:end-7);

%set the variable sized fields to empty
handles.xminm = [];
handles.paramNames = [];
handles.paramsOut = [];
handles.sigma = [];
handles.Cov = [];
handles.pData = [];
handles.pModel = [];
handles.pAllStates = [];
handles.Cov = [];
handles.CovPlot = [];
handles.apprCovFlag = [];
handles.knots = [];
handles.maxDegree = [];
handles.ubR = [];

%initialize obsOx, obsOxCpG and errorsBis errorsOx in case we haven't
    %loaded oxBS and MAB_BS - we need this if bec we cant initiliaze before
    %since we don't know the size of obsMatrix
if (isempty(handles.regionNameOx))
    handles.FileName{2} = '';
    handles.obsData(:,:,2) = zeros(size(handles.obsData(:,:,1)));
    handles.obsDataCpG(:,:,:,2) = zeros(size(handles.obsDataCpG(:,:,:,1)));
    handles.obsData(:,:,3) = zeros(size(handles.obsData(:,:,1)));
end
    
%get the number of days
numOfDays = size(handles.dataPoints, 1);

%if no errors file has been given labels = dataPointsName + 0 | 1 |...
%process = rep, ..., rep and errors are set to 0.005 and 0.07
if isempty(handles.errorsFileName)
    %initialize errors 
    %Bis
    handles.errors = zeros(size(handles.dataPoints, 1), 3, 3);
    handles.errors(:,1,1) = 0.005;
    handles.errors(:,2:3,1) = 0.07;
    handles.errors(:,4,1) = 0.03;
    %oxBis
    handles.errors(:,1,2) = 0.005;
    handles.errors(:,2:3,2) = 0.07;
    handles.errors(:,4,3) = 0.03;
    %mabBis
    handles.errors(:,1,3) = 0.005;
    handles.errors(:,2:3,3) = 0.07;
    handles.errors(:,4,3) = 0.03;
    %initialize labels taking only the first letter (day ->d hour ->h)
    handles.labels = cellstr(repmat(handles.dataPointsName, numOfDays, 1));
    for i=1:numOfDays
        handles.labels{i} = strcat(handles.labels{i}(1), num2str(handles.dataPoints(i))); 
    end
    %initialize process
    a = {'rep'};
    handles.process = repmat(a, numOfDays, 1);
%else replace the nan entries with default values    
else
    %BS, oxBS, mabBS
%     handles.errors(isnan(handles.errors(:,1,:)), 1, :) = 0.005;
%     handles.errors(isnan(handles.errors(:,2:3,:)), 2:3, :) = 0.07;
%     handles.errors(isnan(handles.errors(:,4,:)), 2:3, :) = 0.03;
end

%-------RUN THE ESTIMATION PROCEDURE HERE-------
if get(handles.popupmenu_eff, 'Value') == 1

    [handles.xminm, handles.paramNames, handles.paramsOut, handles.sigma, handles.fmin, handles.pData, handles.pModel, ...
        handles.pAllStates, handles.Cov, handles.CovPlot, handles.apprCovFlag, handles.knots, handles.maxDegree, handles.ubR, handles.model]...
        = UIestimateDSHydroxy(handles.dataPoints, handles.obsData, handles.errors, handles.process, handles.modelFileName);    
else
    [handles.pData, handles.pModel, handles.pAllStates] = UIestimateLevels(handles.dataPoints, handles.obsData, handles.errors);        
end

% handles.errorsMSI,
% handles.obsMSI,
numOfCpGs = size(handles.CpGs, 2);

%set the variable sized fields to empty
handles.xminmCpG = [];
handles.paramsOutCpG = [];
handles.sigmaCpG = [];
handles.CovCpG = [];
handles.pDataCpG = [];
handles.pModelCpG = [];
handles.pAllStatesCpG = [];
handles.CovCpG = [];
handles.CovPlotCpG = [];
handles.apprCovFlagCpG = [];
handles.maxDegreeCpG = [];
handles.knots = [];

for k=1:numOfCpGs
    s = size(handles.obsDataCpG(:,:,k,:));
    obsDataCpGk = reshape(handles.obsDataCpG(:,:,k,:), [s(1), s(2), s(4), s(3)]);
    
    
    if handles.CpGflags(handles.CpGs(k)) == 1
        if get(handles.popupmenu_eff, 'Value') == 1

         [handles.xminmCpG(k,:), ~, handles.paramsOutCpG(k,:), handles.sigmaCpG(k,:), ~, ...
             handles.pDataCpG(:,:,k,:), handles.pModelCpG(:,:,k,:), ...
             handles.pAllStatesCpG(:,:,k), handles.CovCpG(:,:,k), handles.CovPlotCpG(:,:,k), ...
             handles.apprCovFlagCpG(k)]...
             = UIestimateDSHydroxy(handles.dataPoints, obsDataCpGk, handles.errors, handles.process, ... 
                handles.modelFileName);
               
        else          
            [handles.pDataCpG(:,:,k,:), handles.pModelCpG(:,:,k,:), handles.pAllStatesCpG(:,:,k)]...
                = UIestimateLevels(handles.dataPoints, obsDataCpGk, handles.errors);
                      
        end    
    end  
end

%in case of estimating only levels choice fill out the other returned vars 
if get(handles.popupmenu_eff, 'Value') == 2
    handles.xminmCpG = zeros(numOfCpGs, 0);
    handles.paramsOutCpG = zeros(numOfCpGs, 0);
    handles.sigmaCpG = zeros(numOfCpGs, 0);
    handles.CovCpG = zeros(0, 0, numOfCpGs);
    handles.CovPlotCpG = zeros(0, 0, numOfCpGs);
    handles.apprCovFlagCpG = zeros(numOfCpGs, 1);
    
end


%delete msgbox
delete(m);

%update handles
guidata(hObject, handles);

%call UIhydroxy window to plot the results
UIhydroxy('UIdialogBox', handles.figure1);


%% --- Executes on button press in undo_pushbutton6.
function undo_pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to undo_pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dataType = '';
handles.FileName{2} = '';
handles.regionNameOx = [];
handles.obsData(:,:,2) = zeros(size(handles.obsData(:,:,1)));
handles.obsDataCpG(:,:,:,2) = zeros(size(handles.obsDataCpG(:,:,:,1)));

if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
%     set(handles.exportpdf_pushbutton4, 'Enable', 'off');
%     set(handles.exporttxt_pushbutton5, 'Enable', 'off');
end

%clear textbox and tick
set(handles.oxBS_text2, 'Visible', 'off');
cla(handles.tick_axes3);
%set(handles.tick_axes2, 'Visible', 'off');


%clear the axes
% cla(handles.tick_axes2, 'reset');
% cla(handles.tick_axes3, 'reset'); 
% cla(handles.axes3, 'reset'); 
% cla(handles.tick_axes1, 'reset'); 
% cla(handles.axes5, 'reset');
%clear legends
% legend(handles.tick_axes2,'hide');
% legend(handles.tick_axes3,'hide');
% legend(handles.axes3,'hide');
% legend(handles.tick_axes1,'hide');
% legend(handles.axes5,'hide');


guidata(hObject, handles);


% --- Executes on button press in undo_pushbutton7.
function undo_pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to undo_pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.dataType = '';
handles.FileName{3} = '';
handles.regionNameMab = [];
handles.obsData(:,:,3) = zeros(size(handles.obsData(:,:,1)));
handles.obsDataCpG(:,:,:,3) = zeros(size(handles.obsDataCpG(:,:,:,1)));

if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
%     set(handles.exportpdf_pushbutton4, 'Enable', 'off');
%     set(handles.exporttxt_pushbutton5, 'Enable', 'off');
end

%clear textbox and tick
set(handles.mabBS_text, 'Visible', 'off');
cla(handles.tick_axes4);

guidata(hObject, handles);


%% --- Executes on button press in exampleBS_pushbutton7.
function exampleBS_pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to exampleBS_pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = computer;
pathBis = which('/IAP_BS.txt');
fileattrib(pathBis, '-w', 'a'); 

if (strfind(str, 'MAC'))
    %make file only read and open it with TextEdit
    fileattrib(pathBis, '-w', 'a'); 
    system(['open -a TextEdit ' pathBis]);
elseif (strfind(str, 'WIN'))
    %make file only read
    fileattrib(pathBis, '-w', ''); 
    %find the right command for windows
    winopen(pathBis);
    %system(['notepad %s', fullfile(pwd,'L1mdTBis.txt')]);
else 
    %make file only read
    fileattrib(pathBis, '-w', 'a');
    %find the right texteditor for linux
    [~, deskMnger] = system('echo "$XDG_DATA_DIRS" | grep -Eo "xfce|kde|gnome"');
    deskMnger = deskMnger(1:end-1);
    if (strcmp(deskMnger, 'xfce'))
        %find the right command for linux
        system(['mousepad ', pathBis, '&']);
    elseif (strcmp(deskMnger, 'gnome'))
        system(['gedit ', pathBis, '&']);
    elseif (strcmp(deskMnger, 'kde'))
        system(['dolphin ', pathBis, '&']);
    else
        msg = 'No txt editor';
        msgbox('Can not find a text editor');
        error(msg);
    end
end    
    

%% --- Executes on button press in exampleOxBS_pushbutton5.
function exampleOxBS_pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to exporttxt_pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = computer;
pathOx = which('/IAP_oxBS.txt');

if (strfind(str, 'MAC'))
    %make file only read and open it with TextEdit
    fileattrib(pathOx, '-w', 'a'); 
    system(['open -a TextEdit ' pathOx]);
elseif (strfind(str, 'WIN'))
    %make file only readable and open it
    fileattrib(pathOx, '-w', '');
    winopen(pathOx);
    %system(['notepad %s', fullfile(pwd,'L1mdTOx.txt')]);
else
    %make file only read and open it
    fileattrib(pathOx, '-w', 'a');
    [~, deskMnger] = system('echo "$XDG_DATA_DIRS" | grep -Eo "xfce|kde|gnome"');
    deskMnger = deskMnger(1:end-1);
    if (strcmp(deskMnger, 'xfce'))
        %find the right command for linux
        system(['mousepad ', pathOx, '&']);
    elseif (strcmp(deskMnger, 'gnome'))
        system(['gedit ', pathOx, '&']);
    elseif (strcmp(deskMnger, 'kde'))
        system(['dolphin ', pathOx, '&']);
    else
        msg = 'No txt editor';
        msgbox('Can not find a text editor');
        error(msg);
    end   
end    


%% --- Executes on button press in exampleMabBS_pushbutton9.
function exampleMabBS_pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to exampleMabBS_pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = computer;
pathMab = which('/DMR2_mabBS.txt');

if (strfind(str, 'MAC'))
    %make file only read and open it with TextEdit
    fileattrib(pathMab, '-w', 'a'); 
    system(['open -a TextEdit ' pathMab]);
elseif (strfind(str, 'WIN'))
    %make file only readable and open it
    fileattrib(pathMab, '-w', '');
    winopen(pathMab);
    %system(['notepad %s', fullfile(pwd,'L1mdTOx.txt')]);
else
    %make file only read and open it
    fileattrib(pathMab, '-w', 'a');
    [~, deskMnger] = system('echo "$XDG_DATA_DIRS" | grep -Eo "xfce|kde|gnome"');
    deskMnger = deskMnger(1:end-1);
    if (strcmp(deskMnger, 'xfce'))
        %find the right command for linux
        system(['mousepad ', pathMab, '&']);
    elseif (strcmp(deskMnger, 'gnome'))
        system(['gedit ', pathMab, '&']);
    elseif (strcmp(deskMnger, 'kde'))
        system(['dolphin ', pathMab, '&']);
    else
        msg = 'No txt editor';
        msgbox('Can not find a text editor');
        error(msg);
    end   
end    

%% --- Executes on button press in errors_pushbutton9.
function errors_pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to errors_pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    [handles.errorsFileName, handles.PathName] = uigetfile('*.txt', 'Select errors File', ...
    '/Users/kyriakopou/Documents/code/MATLAB/HydroxyMethylation/dataFilesNew/newFormat/' ,'Multiselect', 'off');

if (handles.errorsFileName ~= 0) %if no cancel was pressed
    
    %check the name of errorsfileName (_errors.txt)
    if (size(handles.errorsFileName, 2) < 11 || ~strcmp(handles.errorsFileName(end-10:end), '_errors.txt') )
        msg = 'Not appropriate errors fileName';
        msgbox('Please give appropriate _errors.txt fileName');
        error(msg); 
    else
            
        %set the names of errors file
        file = strcat(handles.PathName, handles.errorsFileName);
        
        %call readErrors 
        [handles.dataPointsName, handles.dataPoints, handles.labels, ...
            handles.errors, handles.process] = UIreadErrors(file);
        
        %initialize errors in case nan is returned for some entries
        defaultErrors = nan(size(handles.errors));
        %Bis
        defaultErrors(:,1,1) = 0.005;
        defaultErrors(:,2:3,1) = 0.07;
        defaultErrors(:,4,1) = 0.03;
        %oxBis
        defaultErrors(:,1,2) = 0.005;
        defaultErrors(:,2:3,2) = 0.07;
        defaultErrors(:,4,2) = 0.03;
        %mabBis
        defaultErrors(:,1,3) = 0.05;
        defaultErrors(:,2:3,3) = 0.07;
        defaultErrors(:,4,3) = 0.03;
        
        %get the indices of nans in the input error matrix and replace them
        %with the default values
        ind = find(isnan(handles.errors));
        handles.errors(ind) = defaultErrors(ind);
        
        %switch on errors_text
        set(handles.errors_text3, 'String', handles.errorsFileName(1:end-4));
        set(handles.errors_text3, 'Visible', 'on');
        
        
        %put the tick image
        axes(handles.tick_axes1);
        Img = imread('tick.png', 'BackgroundColor', [0 0 0]);
        mask = bsxfun(@eq,Img,reshape([0 0 0], 1,1,3));
        handles.tick1 = image(Img,'alphadata', 1-double(all(mask,3)));

        axis off;
        axis image;  
    
        %if BS is already given enable params (in case we change only errors file and we want to run again)
        if ~isempty(handles.regionNameBS)
            set(handles.params_pushbutton3, 'Enable','on');
    %         set(handles.exportpdf_pushbutton4, 'Enable', 'off');
    %         set(handles.exporttxt_pushbutton5, 'Enable', 'off');
        end    
        
    end
    
    
    
    
    %MAYBE ADD SOME CHECKS REGARDING THE NUMBER OF ROWS AND SO ON..
%     %if errors file already given check if #rows of BS agrees with #rows of oxBS
%     %and if not halt and return error msg -> number of data Points differs
%     if (~isempty(handles.FileName{2}))
%         if (size(handles.dataPointsOx, 1) ~= size(handles.dataPoints, 1)) 
%             msg = 'Not appropriate input files';
%             msgbox('Number of dataPoints differs in BS and oxBS files');
%             error(msg); 
%         end
%     else
  
end    

guidata(hObject, handles);

%% --- Executes on button press in errorsEx_pushbutton10.
function errorsEx_pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to errorsEx_pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = computer;
pathErrors = which('/IAP_errors.txt');

if (strfind(str, 'MAC'))
    fileattrib(pathErrors, '-w', 'a'); 
    system(['open -a TextEdit ' pathErrors]);
elseif (strfind(str, 'WIN'))
    fileattrib(pathErrors, '-w', ''); 
    %find the right command for windows
    winopen(pathErrors);
    %system(['notepad %s', fullfile(pwd,'L1mdTBis.txt')]);
else 
    fileattrib(pathErrors, '-w', 'a'); 
    %find the right texteditor for linux
    [~, deskMnger] = system('echo "$XDG_DATA_DIRS" | grep -Eo "xfce|kde|gnome"');
    deskMnger = deskMnger(1:end-1);
    if (strcmp(deskMnger, 'xfce'))
        %find the right command for linux
        system(['mousepad ', pathErrors, '&']);
    elseif (strcmp(deskMnger, 'gnome'))
        system(['gedit ', pathErrors, '&']);
    elseif (strcmp(deskMnger, 'kde'))
        system(['dolphin ', pathErrors, '&']);
    else
        msg = 'No txt editor';
        msgbox('Can not find a text editor');
        error(msg);
    end
end    


%% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
% Hint: get(hObject,'Value') returns toggle state of checkbox1

if get(hObject, 'Value') == 1
    handles.CpGflags(1) = 1;
else
    handles.CpGflags(1) = 0;
    set(handles.selAll, 'value', 0);
end    

 %if BS data already given enable params button
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
numOfCpGs = size(handles.CpGs, 2);
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if get(hObject, 'Value') == 1
    handles.CpGflags(2) = 1;
else
    handles.CpGflags(2) = 0;
    set(handles.selAll, 'value', 0);
end

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end

%if all checkboxes have been checked by hand switch selAll to 1
numOfCpGs = size(handles.CpGs, 2);
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


%% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox3
if get(hObject, 'Value') == 1
    handles.CpGflags(3) = 1;
else
    handles.CpGflags(3) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


%% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox4
if get(hObject, 'Value') == 1
    handles.CpGflags(4) = 1;
else
    handles.CpGflags(4) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


%% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox5
if get(hObject, 'Value') == 1
    handles.CpGflags(5) = 1;
else
    handles.CpGflags(5) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


%% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(6) = 1;
else
    handles.CpGflags(6) = 0; 
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS) 
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(7) = 1;
else
    handles.CpGflags(7) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(8) = 1;
else
    handles.CpGflags(8) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(9) = 1;
else
    handles.CpGflags(9) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

if get(hObject, 'Value') == 1
    handles.CpGflags(10) = 1;
else
    handles.CpGflags(10) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(11) = 1;
else
    handles.CpGflags(11) = 0; 
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(12) = 1;
else
    handles.CpGflags(12) = 0;
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(13) = 1;
else
    handles.CpGflags(13) = 0;
    set(handles.selAll, 'value', 0);
end 

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox6.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if get(hObject, 'Value') == 1
    handles.CpGflags(14) = 1;
else
    handles.CpGflags(14) = 0; 
    set(handles.selAll, 'value', 0);
end    

%if BS data already given enable params button (if it is not)
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable', 'on');
end
%if all checkboxes have been checked by hand switch selAll to 1
if all (handles.CpGflags(handles.CpGs) == 1)
    set(handles.selAll, 'value', 1);
end

guidata(hObject, handles);


% --- Executes on button press in selectAll_checkbox7.
function selectAll_Callback(hObject, eventdata, handles)
% hObject    handle to selectAll_checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%WHY HERE WE HAVE TO GET IT FROM hObject?
handles = guidata(hObject);

% Hint: get(hObject,'Value') returns toggle state of selectAll_checkbox7
if get(hObject,'Value') == 1
    for i = handles.CpGs
%       set(handles.(sprintf('checkbox%d',i)),'value', 1)
        set(handles.checkbox(i),'value', 1)
        handles.CpGflags(i) = 1;
    end    
else
    for i = handles.CpGs
        set(handles.checkbox(i),'value', 0)
        handles.CpGflags(i) = 0;
    end    
end    

%if BS is already given enable params
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
%     set(handles.exportpdf_pushbutton4, 'Enable', 'off');
%     set(handles.exporttxt_pushbutton5, 'Enable', 'off');
end    

guidata(hObject, handles);



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg('Are you sure you want to exit H(O)TA?',...
                 'Close Request Function',...
                 'Yes', 'No', 'Yes');
switch selection,
    case 'Yes',
        %close first UIhydroxy
        hf=findobj('Name', 'H(O)TA Main Window');
        close(hf);
        % destroy the dialogBox now.
        delete(hObject); 
end
 


% --- Executes on selection change in popupmenu_model.
function popupmenu_model_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_model
%if BS is already given enable params
if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
%     set(handles.exportpdf_pushbutton4, 'Enable', 'off');
%     set(handles.exporttxt_pushbutton5, 'Enable', 'off');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_eff.
function popupmenu_eff_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_eff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_eff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_eff
%if BS is already given enable params
% handles = guidata(hObject);

% if get(handles.popupmenu_eff, 'Value') == 1
%     set(handles.popupmenu_model, 'Enable', 'on');
% end
% 
% if get(handles.popupmenu_eff, 'Value') == 2
%     set(handles.popupmenu_model, 'Enable', 'off');
% end

if ~isempty(handles.regionNameBS)
    set(handles.params_pushbutton3, 'Enable','on');
%     set(handles.popupmenu_model, 'Enable', 'on');
end    

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenu_eff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_eff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
