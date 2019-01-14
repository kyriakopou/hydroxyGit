function [dataPointsName, dataPoints, labels, obs, errors, process] = UIreadData(file)
%read Data from the files Bis.txt and Ox.txt  

%define the number of observable states
numObsStates = 4;

%read the observations and the error constants from
%the appropriate .txt file and store in tables
%T = readtable(file, 'Delimiter',',');

try
    T = readtable(file, 'Delimiter',',');
catch
    msg = 'Not appropriate input file';
    msgbox('Please give appropriate input file');
    error(msg) 
end


%control for the input files format
%add checks for the format of each column
if (size(T, 2) ~= 10)
    msg = 'Not appropriate input file';
    msgbox('Please give appropriate input file');
    error(msg); 
    
else
    %vector containing the days of the observations
    dataPoints = table2array(T(:, 1));
    %array containing the labels of the dataPoints
    labels = table2array(T(:, 2));

    %Matrices that contain all the observations
    %taken from bisulfite and Oxbisulfite
    obs = table2array(T(:, 3:numObsStates+2));
    %column vectors of the conversion errors for
    %each time point
    errors = table2array(T(:, numObsStates+3:numObsStates+5));
    %last column defines if the transition btw two dataPoints 
    %is a cell-replication or not
    process = table2array(T(:, numObsStates+6));
    %first variable name (either day or dataPoint)
    dataPointsName = T.Properties.VariableNames(1);
    
    %compute the minimum number of dataPoints of the set
%     minPoints = minDataPoints(obs);
    
    
end

end