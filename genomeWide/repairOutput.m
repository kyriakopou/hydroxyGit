function repairOutput

%read the whole results.sorted.igv file
T_igv = readtable('/Users/kyriakopou/Desktop/results.sorted.txt', 'Delimiter', '\t');   
T_chrm = readtable('/Users/kyriakopou/Desktop/chrmCol.txt', 'Delimiter', '\t');

% T_igv.Start = T_igv.Start +1;
% T_igv.End = T_igv.End -1;

T_igv.Chromosome = T_chrm; 

writetable(T_igv, '/Users/kyriakopou/Desktop/results.sorted.txt', 'Delimiter', '\t');

end