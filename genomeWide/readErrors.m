function [dataPointsName, dataPoints, labels, errors_BS, errors_oxBS, process] = readErrors(file)
%read Data from the files Bis.txt and Ox.txt  

%read the observations and the error constants from
%the appropriate .txt file and store in tables

try
    T = readtable(file, 'Delimiter',',');
catch
    msg = 'Not appropriate input file';
    msgbox('Please give appropriate input file\n See the example file');
    error(msg) 
end


%control for the input files format
%add checks for the format of each column
if (size(T, 2) ~= 9)
    msg = 'Not appropriate input file';
    msgbox('Please give appropriate input file. See the example file');
    error(msg); 
    
else
    %vector containing the days of the observations
    dataPoints = table2array(T(:, 1));
    %array containing the labels of the dataPoints
    labels = table2array(T(:, 2));

    %Matrices that contain all the observations
    %taken from bisulfite and Oxbisulfite
    errors_BS = table2array(T(:, 3:5));
    %column vectors of the conversion errors for
    %each time point
    errors_oxBS = table2array(T(:, 6:8));
    %last column defines if the transition btw two dataPoints 
    %is a cell-replication or not
    process = table2array(T(:, 9));
    %first variable name (either day or dataPoint)
    dataPointsName = T.Properties.VariableNames(1);
    
    %compute the minimum number of dataPoints of the set
%     minPoints = minDataPoints(obs);
    
    
end

end