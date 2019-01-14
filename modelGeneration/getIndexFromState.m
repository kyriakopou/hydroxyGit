function index = getIndexFromState(state, varargin)
%general implementation of encoding a state (string of methylation status
%chars) to a decimal index

numOfSites = length(state);
nVarargs = length(varargin);

if (nVarargs == 0)
    varValues = getVarValues;
elseif (nVarargs == 1)
    varValues = varargin{1};
else
    error('wrong number of input arguments');
end

%order varValues
varValues = fliplr(sort(varValues));
numOfValues = size(varValues, 2);

%check if all characters given in state
for i=1:numOfSites
    if isempty(find(strcmp(varValues, state(i))));     
        error('The given state contains not valid varValues');
    end    
end

indices = zeros(1, numOfSites);
%find the index for the letter on position k (counting from the left)  
%and store it in indices(k) array
for k=1:numOfSites
    indices(k) = find(strcmp(state(k), varValues));
end

%create an enumeration system with base = numOfVars
%and add one to each index to get indices within range [1, numOfVars^(numOfSites-1) + 1]
powersOfNumOfValues = ones(1, numOfSites) * numOfValues;
powers = numOfSites-1:-1:0;
powersOfNumOfValues = powersOfNumOfValues .^ powers;
%index
index = sum((indices-1) .* powersOfNumOfValues) + 1;


end