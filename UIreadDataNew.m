function [dataPointsName, dataPoints, ObsAgg, Obs, CpGs] = UIreadDataNew(file)
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
    msgbox({'Please give appropriate input file'; 'by checking the example file'});
    error(msg) 
end


%control for the input files format
%add checks for the format of each column
if (size(T, 2) ~= 6)
    msg = 'Not appropriate input file';
    msgbox({'Please give appropriate input file'; 'by checking the example file'});
    error(msg); 
    
else
    
    %first variable name (either day or dataPoint)
    dataPointsName = T.Properties.VariableNames{1};
    
    %sort the row wrt to CpG number first
    T = sortrows(T, {'CpG', dataPointsName}, {'ascend', 'ascend'});
    
    %vector containing the days of the observations
    dataPoints = unique(table2array(T(:, 1)));
    %vector containing the CpGs of the region
    CpGs = unique(table2array(T(:, 2)))';
    %number of timePoints
    numOfTimePoints = size(dataPoints, 1);
    %number of CpGs provided
    numOfCpGs = size(CpGs, 2);
    
    obsTemp = table2array(T(:, 3:numObsStates+2));
         
    %go through obsTemp and fill all Obs matrices and ObsAgg
    %WORKS ONLY IF DATA IS ORDERED 
%     for i=1:numOfTimePoints
%         for k=1:numOfCpGs
%             Obs(i,:,k) = obsTemp((i-1)*numOfCpGs + k, :);
%             ObsAgg(i,:) = ObsAgg(i,:) + Obs(i,:,k);
%         end      
%     end
    %Obs is now a 3D matrix (one 2d for each CpG)
    Obs = zeros(numOfTimePoints, numObsStates, numOfCpGs);
    for k=1:numOfCpGs  
        Obs(:,:,k) = obsTemp((k-1)*numOfTimePoints+1:(k-1)*numOfTimePoints+numOfTimePoints, :);
    end    
    ObsAgg = sum(Obs, 3);
    
    
    
    %compute the minimum number of dataPoints of the set
%     minPoints = minDataPoints(obs);
    
    
end

end