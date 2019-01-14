function metropHastingsExample()
%METROPHASTINGS Summary of this function goes here
%   Detailed explanation goes here

rng default  % For reproducibility

delta = .5;
pdf = @(x) normpdf(x);
%pdf should be now the prior * likellihood 

proppdf = @(x,y) unifpdf(y-x,-delta,delta);
proprnd = @(x) x + rand*2*delta - delta;


nsamples = 500;
% x = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd, 'proppdf', proppdf, 'symmetric',1);
x = mhsample(1,nsamples,'pdf',pdf, 'proprnd',proprnd, 'symmetric', 1, 'burnin', 1000);

figure(2);
h = histfit(x);
h(1).FaceColor = [.8 .8 1];

end

