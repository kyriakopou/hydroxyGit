function plotHiddenStatesOverTime(T_out, labels)

%avg hidden states levels over time
mmOverTime = [T_out.(genvarname(strcat('mm_', labels{1}))), T_out.(genvarname(strcat('mm_', labels{2}))), ... 
    T_out.(genvarname(strcat('mm_', labels{3})))]; 

tothOverTime = [T_out.(genvarname(strcat('toth_', labels{1}))), T_out.(genvarname(strcat('toth_', labels{2}))), ... 
    T_out.(genvarname(strcat('toth_', labels{3})))]; 

hemiOverTime = [T_out.(genvarname(strcat('um_', labels{1}))), T_out.(genvarname(strcat('um_', labels{2}))), ... 
    T_out.(genvarname(strcat('um_', labels{3})))]; 

uuOverTime = [T_out.(genvarname(strcat('uu_', labels{1}))), T_out.(genvarname(strcat('uu_', labels{2}))), ... 
    T_out.(genvarname(strcat('uu_', labels{3})))]; 


levelsAll = {mmOverTime; 
            tothOverTime;
            hemiOverTime;
            uuOverTime};

% old colors
% rgbColor = [0.8500  0.3250  0.0980;
%             0.929   0.6940  0.125;
%             0.466   0.674   0.188;
%             0       0.447   0.741];
        
rgbColor = [1.0000    0.4000         0;
            0.9020    0.7059    0.1137;        
            0.466   0.8   0.188;
            0    0.8000    1.0000];
        
%levels over time only for the three days obs        
if all(any(~isnan(levelsAll{1,1})))
    aboxplot(levelsAll, 'Colormap', rgbColor, 'labels', labels);
else
    aboxplot(levelsAll, 'Colormap', rgbColor, 'labels', labels);
end
% legend('mm', 'toth', 'um-mu', 'uu');

end