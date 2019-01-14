function [pinit, samplesBS, samplesOx] = generateRandomData(pTotalHydroxy, EMatrix)

%create multinomial initial distribution for states
%uu um uh uf mu mm mh mf hu hm hh hf fu fm fh ff 
pinit = zeros(1, 16);
%get random probabilities for the distinct hydroxystates that sum up to pHydrox
pVectorHydroxy = rand(1, 5);
pVectorHydroxy = pVectorHydroxy / sum(pVectorHydroxy) * pTotalHydroxy;

%fix hydroxylation and methylation (mm, um, mu levels)
pinit([3, 7, 9:11]) = pVectorHydroxy;
ptotalHydroxy = sum(pVectorHydroxy);

pinit(getIndexFromState('um')) = 0.05;
pinit(getIndexFromState('mu')) = 0.05;
pinit(getIndexFromState('mm')) = 0.1;

ptotalMethyl = sum(pinit(getIndexFromState('um')) + pinit(getIndexFromState('mu')) + pinit(getIndexFromState('mm')));

%remain prob mass to state uu
pinit(getIndexFromState('uu')) = 1 - ptotalMethyl - ptotalHydroxy; 
  
pBis = pinit * EMatrix(:,:,:,1);
pOx = pinit * EMatrix(:,:,:,2);

numOfSamples = 0.5*10^4;
samplesBS = mnrnd(numOfSamples, pBis);
samplesOx = mnrnd(numOfSamples, pOx);


end