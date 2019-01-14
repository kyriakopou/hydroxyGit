function numOfCpGsPerChr = getNumOfCpGsPerChromosome(indFirstCpG, totNumOfCpGs)

numOfChromosomes = size(indFirstCpG, 1);
indLastCpG = indFirstCpG(2:end)-1;
indLastCpG(numOfChromosomes) = totNumOfCpGs;

numOfCpGsPerChr = indLastCpG - indFirstCpG + 1;


end