function plotHiddenStatesOverTime(T_out, labels)

%avg hidden states levels over time
mmOverTime = [T_out.mm_d0, T_out.mm_d3, T_out.mm_d6];
tothOverTime = [T_out.toth_d0, T_out.toth_d3, T_out.toth_d6];
hemiOverTime = [T_out.um_d0, T_out.um_d3, T_out.um_d6];
uuOverTime = [T_out.uu_d0, T_out.uu_d3, T_out.uu_d6];

levelsAll = {mmOverTime; 
            tothOverTime;
            hemiOverTime;
            uuOverTime};

%levels over time
rgbColor = [0.8500  0.3250  0.0980;
            0.929   0.6940  0.125;
            0.466   0.674   0.188;
            0       0.447   0.741];


aboxplot(levelsAll, 'Colormap', rgbColor, 'labels', labels);
% legend('mm', 'toth', 'um-mu', 'uu');

end