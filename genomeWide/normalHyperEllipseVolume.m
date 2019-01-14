function vol = normalHyperEllipseVolume(Sigma, confLevel)
%computes the hypervolume of the hyperellipse of the normal distribution
%with covariance matrix Sigma for a given confLevel of the estimated params

k = size(Sigma, 1);
chi2_crit = chi2inv(confLevel, k);

deter = det(Sigma);
vol = 2*pi^(k/2)/(k*gamma(k/2)) * chi2_crit^(k/2)  * sqrt(deter);


end