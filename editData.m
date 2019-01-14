function editData(dataPath)
%read Data from the files Bis.txt and Ox.txt  

%define the number of observable states
numObsStates = 4;

%read the observations and the error constants from
%the appropriate .txt file and store in tables
file1 = strcat(dataPath, '_BS.txt');
file2 = strcat(dataPath, '_oxBS.txt');

Tbis = readtable(file1, 'Delimiter', ',');
Tox = readtable(file2, 'Delimiter', ',');

%vector containing the days of the observations
dataPoints = table2array(Tbis(:, 1));
%array containing the labels of the dataPoints
label = cell(size(dataPoints, 1), 1);
process = cell(size(dataPoints, 1), 1);

for i=1:size(dataPoints, 1) 
    label{i} = strcat('day', int2str(dataPoints(i)));
    if (i ~= 1)
        process{i} = 'rep';
    else
        process{i} = '-';
    end    
end

label = table(label);
process = table(process);

TmodBis = [Tbis(:,1), label, Tbis(:,2:8), process]; 
TmodOx = [Tox(:,1), label, Tox(:,2:8), process]; 

TmodBis.Properties.VariableNames{7} = 'CError';
TmodBis.Properties.VariableNames{8} = 'mCError';
TmodBis.Properties.VariableNames{9} = 'hmCError';

TmodOx.Properties.VariableNames{7} = 'CError';
TmodOx.Properties.VariableNames{8} = 'mCError';
TmodOx.Properties.VariableNames{9} = 'hmCError';


%write these tables to files
writetable(TmodBis, file1);
writetable(TmodOx, file2);

end