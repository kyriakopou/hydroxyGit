function lambdaCoef(params, Cov, fileID)

%some necessary moments
%second
secondb0Mm_b0Md = Cov(1,3) + params(1) * params(3);
secondb1Mm_b1Md = Cov(2,4) + params(2) * params(4);
secondb0Mm_b1Md = Cov(1,4) + params(1) * params(4);
secondb1Mm_b0Md = Cov(2,3) + params(2) * params(3);

%third
thirdb0Mm2_b0Md = params(1)^2*params(3) + 2*Cov(1,3)*params(1) + Cov(1,1)*params(3);
thirdb0Mm_b0Md2 = params(3)^2*params(1) + 2*Cov(1,3)*params(3) + Cov(3,3)*params(1);



%fourth
fourthb0Mm2_b0Md2 = params(1)^2*params(3)^2 + Cov(1,1)*Cov(3,3) + 2*Cov(1,3)^2 + ...
    4*Cov(1,3)*params(1)*params(3) + Cov(1,1)*params(3)^2 + Cov(3,3)*params(1)^2;

fourthb1Mm2_b1Md2 = params(2)^2*params(4)^2 + Cov(2,2)*Cov(4,4) + 2*Cov(2,4)^2 + ...
    4*Cov(2,4)*params(2)*params(4) + Cov(2,2)*params(4)^2 + Cov(4,4)*params(2)^2;

fourthb0Mm2_b1Md2 = params(1)^2*params(4)^2 + Cov(1,1)*Cov(4,4) + 2*Cov(1,4)^2 + ...
    4*Cov(1,4)*params(1)*params(4) + Cov(1,1)*params(4)^2 + Cov(4,4)*params(1)^2;

fourthb1Mm2_b0Md2 = params(2)^2*params(3)^2 + Cov(2,2)*Cov(3,3) + 2*Cov(2,3)^2 + ...
    4*Cov(2,3)*params(2)*params(3) + Cov(2,2)*params(3)^2 + Cov(3,3)*params(2)^2;



%variances of products refer to the appendix
varb0Mm_b0Md = fourthb0Mm2_b0Md2 - secondb0Mm_b0Md^2;
varb1Mm_b1Md = fourthb1Mm2_b1Md2 - secondb1Mm_b1Md^2;
varb0Mm_b1Md = fourthb0Mm2_b1Md2 - secondb0Mm_b1Md^2;
varb1Mm_b0Md = fourthb1Mm2_b0Md2 - secondb1Mm_b0Md^2;


Covb0Mm2_b0Md = thirdb0Mm2_b0Md - params(1)*secondb0Mm_b0Md;
Covb0Mm_b0Md2 = thirdb0Mm_b0Md2(1,:) - params(3)*secondb0Mm_b0Md;


%Covb1Mm2_b0Md = thirdb1Mm2_b0Md(1,:) - params(2)*secondb1Mm_b0Md;
% Covb1Md2_b0Mm = 
% Covb1Mm_b0Mm_b1Md = 
% Covb1Md_b1Mm_b0Md = 
% Covb1Mm_b0Md_b0Mm_b1Md = 

%the total methylation of hemimethylated sites (lambda)
b0Lambda = params(1) + params(3) - params(1)*params(3);
b1Lambda = params(2) + params(4) - params(1)*params(4) - params(2)*params(3);
b2Lambda = - params(2)*params(4);

%sum all the previous to get the variance of lambda
varb0Lambda = Cov(1,1) + Cov(3,3) + varb0Mm_b0Md  + 2*(Cov(1,3) - Covb0Mm2_b0Md - Covb0Mm_b0Md2);
% varb1Lambda = Cov(2,2) + Cov(4,4) + varb0Mm_b1Md + varb1Mm_b0Md + 2*(Cov(2,4) - Covb1Mm2_b0Md - ...
%             Covb1Mm_b0Mm_b1Md - Covb1Md_b1Mm_b0Md - Covb1Md2_b0Mm + Covb1Mm_b0Md_b0Mm_b1Md);
varb2Lambda = varb1Mm_b1Md;

fprintf(fileID, 'b0_l = %d\n', b0Lambda);
fprintf(fileID, 'b1_l = %d\n', b1Lambda);
fprintf(fileID, 'b2_l = %d\n', b2Lambda);
fprintf(fileID, '\n');


end