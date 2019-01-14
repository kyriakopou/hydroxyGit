function T_final = preprocessData(dataFilesPath, maxObsDaysVector)

numOfObsDays = length(maxObsDaysVector);

%if T_final is already stored 
% inputMatFile = strcat(dataFilesPath,'/T_test.mat');
% inputMatFile = strcat('T_input.mat');
inputMatFile = strcat(dataFilesPath, '/T_input.mat');

if exist(inputMatFile, 'file') == 2
    load(inputMatFile);
else
    %read the genome wide data
    T = cell(numOfObsDays, 1);

    for i=1:numOfObsDays

        %get the name of the files for each dayfor
        fileNameBS = strcat(dataFilesPath, '/d', int2str(maxObsDaysVector(i)), '_BS.txt');
        T_BS = readtable(fileNameBS, 'Delimiter', '\t');
        
        %remove unnecessary columns
        T_BS.strand = [];
        T_BS.total_cov = [];	
        T_BS.meth_cov = [];	
        T_BS.unmeth_cov = [];
        T_BS.loci_type = [];
        
        
        fileNameOx = strcat(dataFilesPath, '/d', int2str(maxObsDaysVector(i)), '_oxBS.txt');
        %if there is hydroxylation file
        if exist(fileNameOx, 'file') == 2

            T_oxBS = readtable(fileNameOx, 'Delimiter', '\t');

            T_oxBS.strand = [];	
            T_oxBS.total_cov = [];
            T_oxBS.meth_cov = [];	
            T_oxBS.unmeth_cov = [];	
            T_oxBS.loci_type = [];

            %get all the common rows between BS and oxBS for each day
            %i.e., for a day we have BS data we should also have oxBS 
            T{i} = innerjoin(T_BS, T_oxBS, 'Keys', {'chr', 'position'});

        else
            T{i} = T_BS; 
        end
    end
    %also consider rows with less than three days
    T_temp = outerjoin(T{1}, T{2}, 'Keys', {'chr', 'position'}, 'MergeKeys', true);
    T_final = outerjoin(T_temp, T{3}, 'Keys', {'chr', 'position'}, 'MergeKeys', true);
    
    %remove the days with total number of obervations per day and
    %experiments less than 5 samples
    T_final = filteringTable(T_final);
        
    %store T_final in .mat file for next time
    save(inputMatFile, 'T_final');
    
    %delete unnecessary data structures
    clear T_temp T_BS T_oxBS;
    for i=1:numOfObsDays
        clear T{i};
    end
    
end    

% T_final = T_test;
    
%lets consider only rows with observations at all three days
%to get the index vectors
% T_OneAndTwo = innerjoin(T{1}, T{2}, 'Keys', {'chr', 'position'});
% T_final = innerjoin(T_OneAndTwo, T{3}, 'Keys', {'chr', 'position'});


%sort the rows of T_final
% T_final = sortrows(T_final, {'chr', 'position'}, {'ascend', 'ascend'});

%if T_final.chr column is not string make it
if ~iscellstr(T_final.chr)
    T_final.chr = num2str(T_final.chr);
end    


end