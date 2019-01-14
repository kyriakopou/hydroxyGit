function KL = kullbackLeibler(P, Q)
%computes the total kullback leibler distance between the data
%distribution matrix P (distribution over many time points)
%and the model prediction distribution matrix Q

KL = nansum(nansum(P.*log(P./Q), 2));
% KL = nanmean(KL(~isinf(KL)));


end