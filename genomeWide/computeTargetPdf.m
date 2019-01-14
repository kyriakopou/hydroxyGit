function logPdf = computeTargetPdf(x, priorPdf, logLikelihoodPdf)
%Computes the logPdf of the target distribution (proportional to the posterior pdf)

% compute the density of the given prior and likelihood
logPrior = log(priorPdf(x));
logLik = logLikelihoodPdf(x);

%multiply prior with likelihood (or add the logs)
logPdf = logPrior + logLik;

%print warning sign
% if logPdf == 0 || isnan(logPdf) || isinf(logPdf)
%     disp('stop');
% end
    
end