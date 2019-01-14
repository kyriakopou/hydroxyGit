function createNewInputFiles(regionName, CpGs)
%CREATENEWINPUTFILES Summary of this function goes here
%   Detailed explanation goes here

%Create first the error files
%read the observations and the error constants from
%the appropriate .txt file and store in tables
file1 = strcat('Wt/', regionName, '_BS.txt');
file2 = strcat('Wt/', regionName, '_oxBS.txt');
Tbis = readtable(file1, 'Delimiter',',');
Tox = readtable(file2, 'Delimiter',',');

%change errorNames
Tbis.Properties.VariableNames{'CError'} = 'CErr_BS';
Tbis.Properties.VariableNames{'mCError'} = 'mCErr_BS';
Tbis.Properties.VariableNames{'hmCError'} = 'hmCErr_BS';
Tox.Properties.VariableNames{'CError'} = 'CErr_oxBS';
Tox.Properties.VariableNames{'mCError'} = 'mCErr_oxBS';
Tox.Properties.VariableNames{'hmCError'} = 'hmCErr_oxBS';

TErrors = [Tbis(:,1:2), Tbis(:,7:9), Tox(:,7:10)];
numOfDays = size(TErrors, 1);

%write to the errors file
% fileID = fopen(strcat('newFormat/', regionName, '_errors.txt'), 'wt');
writetable(TErrors, strcat('newFormat/', regionName, '_errors'), 'FileType', 'text');
% fclose(fileID);

T_BS = cell2table(cell(0,6), 'VariableNames', {'day', 'CpG', 'TT', 'TC', 'CT', 'CC'});
T_oxBS = cell2table(cell(0,6), 'VariableNames', {'day', 'CpG', 'TT', 'TC', 'CT', 'CC'});

for i = CpGs
    file1 = strcat('singleCpGs/', regionName, 'CpG', num2str(i), '_BS.txt');
    file2 = strcat('singleCpGs/', regionName, 'CpG', num2str(i), '_oxBS.txt');
    %BS oxBS tables for CpG(i)
    Tbis = readtable(file1);
    Tox = readtable(file2);
    
    CpG = repmat(i, numOfDays, 1);
    TCpGCol = table(CpG);
    
    T_BS = [T_BS; Tbis(:,1), TCpGCol, Tbis(:,3:6)];
    T_oxBS = [T_oxBS; Tox(:,1), TCpGCol, Tox(:,3:6)];
    
end

%order the tables wrt to days and then CpG numbers
T_BS = sortrows(T_BS, {'day', 'CpG'}, {'ascend', 'ascend'});
T_oxBS = sortrows(T_oxBS, {'day', 'CpG'}, {'ascend', 'ascend'});

%write to corresponding files 
writetable(T_BS, strcat('newFormat/', regionName, '_BS'), 'FileType', 'text');

writetable(T_oxBS, strcat('newFormat/', regionName, '_oxBS'), 'FileType', 'text');

end

