function matrixFileName = produceQMatrixFile(maxDegree, knots)

numOfHiddenStates = 16;

numOfRates = 5;
numOfKnots = size(knots, 2);

%Check if fileName exists. Otherwise produce the file
a = strsplit(num2str(knots));
matrixFileName = strcat('QMatrix_deg_', int2str(maxDegree));
if ~isempty(knots)
    matrixFileName = strcat(matrixFileName, '_', a{:});
end    
fileName = strcat('~/Documents/code/MATLAB/Hydroxymethylation/modelGeneration/matrices/process/active/', ...
    matrixFileName, '.m');

if ~exist(fileName, 'file')
    fprintf('Producing necessary QmatrixFile. It might take a while... \t');
    
    numOfCoef = maxDegree + 1 + numOfKnots;
    numOfParams = numOfRates * numOfCoef + 1;

    Q = sym(zeros(numOfHiddenStates, numOfHiddenStates));
    derQ = sym(zeros(numOfHiddenStates, numOfHiddenStates, numOfParams));
    secDerQ = sym(zeros(numOfHiddenStates, numOfHiddenStates, numOfParams, numOfParams));

    syms m_m m_d h f d t

    %DEFINE THE TRANSITIONS OF THE CTMC
    Q(1,[2,5]) = m_d;
    Q(1,1) = -2*m_d;
    Q(2,3) = h;
    Q(2,6) = m_d;
    Q(2,2) = -h-m_d;
    Q(3,4) = f;
    Q(3,7) = m_d;
    Q(3,3) = -m_d -f;
    Q(4,1) = d;
    Q(4,8) = m_d;
    Q(4,4) = -d-m_d;
    Q(5,6) = m_d;
    Q(5,9) = h;
    Q(5,5) = -h-m_d;
    Q(6,7) = h;
    Q(6,10) = h;
    Q(6,6) = -2*h;
    Q(7,8) = f;
    Q(7,11) = h;
    Q(7,7) = -h-f;
    Q(8,5) = d;
    Q(8,12) = h;
    Q(8,8) = -h-d;
    Q(9,13) = f;
    Q(9,10) = m_d;
    Q(9,9) = -f-m_d;
    Q(10,14) = f;
    Q(10, 11) = h;
    Q(10,10) = -h-f;
    Q(11,12) = f;
    Q(11,15) = f;
    Q(11,11) = -2*f;
    Q(12,16) = f;
    Q(12,9) = d;
    Q(12,12) = -d-f;
    Q(13,1) = d;
    Q(13,14) = m_d;
    Q(13,13) = -d-m_d;
    Q(14,2) = d;
    Q(14,15) = h;
    Q(14,14) = -d-h;
    Q(15,3) = d;
    Q(15,16) = f;
    Q(15,15) = -d-f;
    Q(16,4) = d;
    Q(16,13) = d;
    Q(16,16) = -2*d;


    %first derivatives (d/db_0 for every efficiency)
    derQ(:,:,1) = diff(Q, m_m);
    derQ(:,:,1+numOfCoef) = diff(Q, m_d);
    derQ(:,:,1+2*numOfCoef) = diff(Q, h);
    derQ(:,:,1+3*numOfCoef) = diff(Q, f);
    derQ(:,:,1+4*numOfCoef) = diff(Q, d);

    %sec derivatives (d^2/db_0^2 for every efficiency)
    secDerQ(:,:,1,1) = diff(Q, m_m, 2);
    secDerQ(:,:,1+numOfCoef) = diff(Q, m_d, 2);
    secDerQ(:,:,1+2*numOfCoef, 1+2*numOfCoef) = diff(Q, h, 2);
    secDerQ(:,:,1+3*numOfCoef, 1+3*numOfCoef) = diff(Q, f, 2);
    secDerQ(:,:,1+4*numOfCoef, 1+4*numOfCoef) = diff(Q, d, 2);

    polyTerms = t.^(0:maxDegree);
    truncTerms = heaviside(t-knots) .* (t-knots).^(maxDegree);
    baseTerms = [polyTerms, truncTerms];
    baseTermsSecDer = reshape(baseTerms, [numOfCoef 1]) * baseTerms;

    for i=1:numOfRates

        rateInd = (i-1)*numOfCoef;
        %---FIRST DERIVATIVE-----
        for j1=1:numOfCoef        
            derQ(:,:,rateInd+j1) = baseTerms(j1) .* derQ(:,:,rateInd+1);
            %----SECDERIVATIVE-----
            for j2=1:numOfCoef
                secDerQ(:,:,rateInd+j1,rateInd+j2) =  baseTermsSecDer(j1, j2) * secDerQ(:,:,rateInd+1,rateInd+1);
            end
        end

    end    

    %
    matlabFunction(Q, derQ, secDerQ, 'File', fileName, 'Outputs', {'M', 'derM', 'secDerM'});
    disp('done!')
end