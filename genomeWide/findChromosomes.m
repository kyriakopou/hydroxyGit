function [uniqueChromosomes, indFirstCpG] = findChromosomes(T)

[uniqueChromosomes, indFirstCpG, ~] = unique(T(:,1));
uniqueChromosomes = table2cell(uniqueChromosomes);


end