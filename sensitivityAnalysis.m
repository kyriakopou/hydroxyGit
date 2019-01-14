function sensitivityAnalysis(xminm, p0, dataPoints, obsBis, obsOx, EBis, EOx, process, pAllStates, paramNames, regionName)

Plus1PerCentMethyl = zeros(size(paramNames, 2), 1);
Minus1PerCentMethyl = zeros(size(paramNames, 2), 1);
Plus1PerCentHydroxy = zeros(size(paramNames, 2), 1);
Minus1PerCentHydroxy = zeros(size(paramNames, 2), 1);

%methylation states
um = sum(pAllStates(:,2:3), 2);
mm = pAllStates(:,8);
totMethylLevel = um+mm;
%hydroxylation states
uh = sum(pAllStates(:,4:5), 2);
hm = sum(pAllStates(:,6:7), 2);
hh = sum(pAllStates(:,9), 2);
totHydroxyLevel = uh + hm + hh;

k = 1;
for param = xminm   
    for i=1:2
        if (i == 1)
            xminm(k) = param * 1.01;
        else
            xminm(k) = param * 0.99;
        end    
            [~, ~, pBis, pOx, pBisModel, pOxModel, pAllStatesPert] = DSHydroxyEmbryo(xminm, p0, dataPoints, obsBis, obsOx, EBis, EOx, process);
            umPert = sum(pAllStatesPert(:,2:3), 2);
            mmPert = pAllStatesPert(:,8);
            
            uhPert = sum(pAllStatesPert(:,4:5), 2);
            hmPert = sum(pAllStatesPert(:,6:7), 2);
            hhPert = pAllStatesPert(:,9);
            
            totMethylLevelPert = umPert + mmPert;
            totHydroxyLevelPert = uhPert + hmPert + hhPert;
            maxDiffHydroxyLevel = max(abs(totHydroxyLevelPert - totHydroxyLevel));
            maxDiffMethylLevel = max(abs(totMethylLevelPert - totMethylLevel));
            
        if (i == 1)
            Plus1PerCentMethyl(k) = maxDiffMethylLevel;
            Plus1PerCentHydroxy(k) = maxDiffHydroxyLevel;
        else
            Minus1PerCentMethyl(k) = maxDiffMethylLevel;
            Minus1PerCentHydroxy(k) = maxDiffHydroxyLevel;        
        end        
                    
    end
        
    k = k + 1;
end

%throw out the demethylation params
Plus1PerCentMethyl = [Plus1PerCentMethyl(1:6); Plus1PerCentMethyl(9)];
Minus1PerCentMethyl = [Minus1PerCentMethyl(1:6); Minus1PerCentMethyl(9)];
Plus1PerCentHydroxy = [Plus1PerCentHydroxy(1:6); Plus1PerCentHydroxy(9)];
Minus1PerCentHydroxy = [Minus1PerCentHydroxy(1:6); Minus1PerCentHydroxy(9)];


T = table(Plus1PerCentMethyl, Minus1PerCentMethyl, Plus1PerCentHydroxy, Minus1PerCentHydroxy, 'RowNames', paramNames);

fileName = strcat('/Users/kyriakopou/Desktop/plosRevision/sensitivity/', regionName, 'Sens', '.txt');
writetable(T, fileName, 'delimiter', ' ', 'WriteRowNames', 1);

% disp(T);


end