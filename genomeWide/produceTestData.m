function T_test =  produceTestData(inputMatFile)

if exist(inputMatFile, 'file') == 2
    load(inputMatFile);
else
    T_final = cell(3, 26);
end

%if it is not string make it
if ~iscellstr(T_final.chr)
    T_final.chr = num2str(T_final.chr);
end    
T_final.chr = strcat('chr', T_final.chr);
%split the data in different chromosomes
[uniqueChromosomes, indFirstCpG, ~] = unique(T_final(:,1));
uniqueChromosomes = table2cell(uniqueChromosomes);
numOfChromosomes = size(uniqueChromosomes, 1);
totNumOfCpGs = size(T_final, 1);

T_chr = cell(numOfChromosomes, 1);

for k=1:numOfChromosomes
    
    if k < numOfChromosomes
        T_chr{k} = T_final(indFirstCpG(k):indFirstCpG(k+1)-1, :);
    else
        T_chr{k} = T_final(indFirstCpG(k):totNumOfCpGs, :);
    end
    
    T_chr{k} = T_chr{k}(1:min(2000,size(T_chr{k})), :);
    
end

T_test = vertcat(T_chr{:});

save 

end