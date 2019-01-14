function [CV_KL, CV_BC, CV_LIK] = LOOCV(dataPath, dataType, assumption, CpGs)
%PERFORMS LEAVE ONE OUT CROSS VALIDATION FOR A CERTAIN REGION (given in dataPath)

CV_KL = 0;
CV_BC = 0;
CV_LIK = 0;
for i=1:size(CpGs,2)
    
    testCpG = CpGs(i);
    trainCpG = [CpGs(1:i-1), CpGs(i+1:end)];
    
    [~, ~, LIK, KL, BC] = estimateDSHydroxy_LOOCV(dataPath, dataType, assumption, trainCpG, testCpG);
    
    CV_LIK = CV_LIK + LIK;
    CV_KL = CV_KL + KL;
    CV_BC = CV_BC + BC;
    
end
    CV_KL = CV_KL / (size(CpGs,2));
    CV_BC = CV_BC / (size(CpGs,2));
    CV_LIK = CV_LIK / (size(CpGs,2));
    
    fprintf('The average log-likelihood is %d \n', CV_LIK);
    fprintf('The average KL divergence is %d \n', CV_KL);
    fprintf('The average BC distance is %d \n', CV_BC);
   

end