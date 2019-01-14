function [dataPointsName, dataPoints, labels, errors, process] = UIreadErrors(file)
%read Data from the files Bis.txt and Ox.txt  

%read the observations and the error constants from
%the appropriate .txt file and store in tables
%T = readtable(file, 'Delimiter',',');

try
    T = readtable(file, 'Delimiter',',');
catch
    msg = 'Not appropriate input file';
    msgbox('Please give appropriate input file\n See the example file');
    error(msg) 
end


%control for the input files format
%add checks for the format of each column
if (size(T, 2) ~= 9 && size(T, 2) ~= 16)
    msg = 'Not appropriate input file';
    msgbox('Please give appropriate input file. See the example file');
    error(msg);   
else
    %vector containing the days of the observations
    dataPoints = table2array(T(:, 1));
    %numOfDataPoints equals the number of rows of the file
    numOfDataPoints = size(dataPoints, 1);
    %array containing the labels of the dataPoints
    labels = table2array(T(:, 2));
    %first variable name (either day or dataPoint)
    dataPointsName = T.Properties.VariableNames(1);
    %initialize errors array
%     errors = nan(numOfDataPoints, 15);
    errors = nan(numOfDataPoints, 4, 3);
    %THIS NOW WRITTEN BY HAND TO CONSIDER AS VALID THE errorFiles that
    %contain the errors of either two or three experiments
    if (size(T, 2) == 16)
        
        %Matrices that contain all the observations
        %taken from bisulfite and Oxbisulfite
        errors(:,:,1) = table2array(T(:, 3:6));
        %column vectors of the conversion errors for
        %each time point
        errors(:,:,2) = table2array(T(:, 7:10));
        %last entry should be MAB_error
        mabErrors = table2array(T(:, 11:15));
        errors(:,1,3) = mabErrors(:,5).*(1-mabErrors(:,1)) + (1-mabErrors(:,5)).*mabErrors(:,2);
        errors(:,2:4,3) = mabErrors(:,2:4);

        %last column defines if the transition btw two dataPoints 
        %is a cell-replication or not
        process = table2array(T(:, 16));


        %compute the minimum number of dataPoints of the set
    %     minPoints = minDataPoints(obs);
    else
        
        %Matrices that contain all the observations
        %taken from bisulfite and Oxbisulfite
        errors(:,1:3,1) = table2array(T(:, 3:5));
        %column vectors of the conversion errors for
        %each time point
        errors(:,1:3,2) = table2array(T(:, 6:8));
        
        %last column defines if the transition btw two dataPoints 
        %is a cell-replication or not
        process = table2array(T(:, 9));

    end
%     end     
end

end