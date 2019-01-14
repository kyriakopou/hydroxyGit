function testMLEPrecision()
%it tests the precision of out 5hmC levels and compares to the subtraction
%method

%fix errors
%Bis
errors(:,1,1) = 0.005;
errors(:,2:3,1) = 0.07;
errors(:,4,1) = 0.03;
%oxBis
errors(:,1,2) = 0.005;
errors(:,2:3,2) = 0.07;
errors(:,4,2) = 0.03;

%get error matrices
EMatrix(:,:,:,1) = myConvErrorsBisAct(1-errors(1,:,1));
EMatrix(:,:,:,2) = myConvErrorsOxAct(1-errors(1,:,2));

%hydroxy levels to be tested
totalHydroxy = 0:0.025:0.5;
%num of repeats per level
numOfRepeats = 10^2;

hydroxyAbsErrorTotal = zeros(length(totalHydroxy), numOfRepeats);
hydroxyAbsErrorOneStrand = zeros(length(totalHydroxy), numOfRepeats);
hydroxyRelError = zeros(length(totalHydroxy), numOfRepeats);
hydroxyAbsSubtr = zeros(length(totalHydroxy), numOfRepeats);

parfor k = 1:length(totalHydroxy)
    disp(k);
    for i = 1:numOfRepeats
        
        [pinit, samplesBS, samplesOx] = generateRandomData(totalHydroxy(k), EMatrix);
        pInitHydroxyStates = pinit([3, 7, 9:11]);
%         pTotalHydroxyInit = sum(pHydroxyInit);

        samples = cat(3, samplesBS, samplesOx);
        
        %estimate the initial distribution using MLE
        %compute initial probability distribution p0
        %reshape s.t. the third dimension goes over experiments
        % EMatrix = cat(4, EBis, EOx, EMAB);
        E = EMatrix(:,:,1,:);
        EMatrixForNullTime = reshape(E, [size(E, 1), size(E,2), size(E, 4), size(E, 3)]);
        [pMLE, ~, ~] = maxLikelihood(samples, EMatrixForNullTime);
        pMLEHydroxyStates = pMLE([3, 7, 9:11]);
%         pMLETotalHydroxy = sum(pMLEHydroxy);
        
        %probability to get an h out of p0
        prob_5hm_MLE = sum(pMLEHydroxyStates(1:end-1)/2) + pMLEHydroxyStates(end);
        prob_5hmc_init = sum(pInitHydroxyStates(1:end-1)/2) + pInitHydroxyStates(end);

%         hydroxyAbsErrorTotal(k, i) = abs(pTotalHydroxyInit - pMLETotalHydroxy);
        hydroxyAbsErrorOneStrand(k, i) = abs(prob_5hmc_init - prob_5hm_MLE);

%         hydroxyRelError(k, i) = hydroxyAbsError(k, i) / pTotalHydroxyInit;

        %estimate the error doing simple subtraction of frequencies
        bsCytosinesFreq = (2*samplesBS(4) + sum(samplesBS(2:3))) / (2*sum(samplesBS));
        oxCytosinesFreq = (2*samplesOx(4) + sum(samplesOx(2:3))) / (2*sum(samplesOx));
        hydroxySubtractEst = bsCytosinesFreq - oxCytosinesFreq;  
        
        hydroxyAbsSubtr(k, i) = abs(hydroxySubtractEst - prob_5hmc_init);
    end
    
end

fig = figure();
plot(totalHydroxy, [mean(hydroxyAbsErrorOneStrand, 2), mean(hydroxyAbsSubtr, 2)], ...
    'LineWidth', 1, 'LineStyle', '--', 'Marker', 'x');

legend('mle', 'subtraction');
xlabel('5hmC level')
ylabel('5hmC estimate mean abs. error');

saveas(fig, '~/Documents/code/scripts/GW_HOTA/plot.pdf')

end
