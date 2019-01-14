function plotAvgEfficiencies(T_out, labels)

rgbColor = [0.8500    0.3250    0.0980;
            0         0.4470    0.7410;
            0.9290    0.6940    0.1250];

      
effAll = {[T_out.(genvarname(strcat('maint_', labels{1}))), T_out.(genvarname(strcat('maint_', labels{2}))), ...
    T_out.(genvarname(strcat('maint_', labels{3})))]; ...
    [T_out.(genvarname(strcat('deNovo_', labels{1}))), T_out.(genvarname(strcat('deNovo_', labels{2}))), ...
    T_out.(genvarname(strcat('deNovo_', labels{3})))]; ...
    [T_out.(genvarname(strcat('hydroxy_', labels{1}))), T_out.(genvarname(strcat('hydroxy_', labels{2}))), ...
    T_out.(genvarname(strcat('hydroxy_', labels{3})))]};
    
% if any(~isnan(effAll{1,1}))
    aboxplot(effAll, 'Colormap', rgbColor, 'labels', labels);
% end
% legend('\mu_m', '\mu_d', '\eta');
ylim([0,1]);

end