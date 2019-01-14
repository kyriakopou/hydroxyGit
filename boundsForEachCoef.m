function [ubVec, lbVec] = boundsForEachCoef(ubR, lastDay, maxDegree, knots)
%provide the bounds for all coefficients of an efficiency function given its upper bound 

polyTermsInUb = ubR ./ lastDay.^(0:maxDegree);
truncTermsInUb = ubR ./ (lastDay-knots).^maxDegree; 
ubVec = [polyTermsInUb, truncTermsInUb];
lbVec = [0, -ubVec(2:end)];

end