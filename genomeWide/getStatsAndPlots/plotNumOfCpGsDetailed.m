function plotNumOfCpGsDetailed(numOfCpGsPerNumOfTimePoints)
%PLOTNUMOFCPGSDETAILED Summary of this function goes here
%   Detailed explanation goes here

numOfDays = size(numOfCpGsPerNumOfTimePoints, 2);
b = bar(numOfCpGsPerNumOfTimePoints);
% xlim([0.25, 1.75]);
colorMatrix = bone(numOfDays);
% for i=1:numOfDays
%     b(i).FaceColor = colorMatrix(i, :);
% end
b.FaceColor = 'flat';
b.CData(1,:) = [0.3010    0.7450    0.9330];
b.CData(2,:) = [0.3010    0.7450    0.9330];
b.CData(3,:) = [.5 0 .5];
% ylabel('numOfCpGs', 'FontSize', 10);
 
end

