function printSamplesPerChrTable(T_in)

%avg num of observations per chromosome and day
inputVarsBS = {'bs0', 'bs3', 'bs6'};
inputVarsOx = {'ox0', 'ox3', 'ox6'};

T_meanChr_bs = varfun(@nanmean, T_in, 'InputVariables', inputVarsBS,...
    'GroupingVariables', 'chr');

T_meanChr_ox = varfun(@nanmean, T_in, 'InputVariables', inputVarsOx,...
    'GroupingVariables', 'chr');

disp(T_meanChr_bs);
disp(T_meanChr_ox);


end