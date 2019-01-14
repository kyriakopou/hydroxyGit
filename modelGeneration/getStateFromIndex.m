function state = getStateFromIndex(index, varargin)
%return the index of a state. Two optional arguments can be given,
% namely the variableSet and the numOfSites. 


nVarargs = length(varargin);
if nVarargs == 0
    varValues = getVarValues;
    numOfSites = getNumOfSites;
elseif nVarargs == 1
    varValues = varargin{1};
    numOfSites = getNumOfSites;
elseif  nVarargs == 2
    varValues = varargin{1};
    numOfSites = varargin{2};
else 
    error('Not the correct number of arguments')
end

%order the variable's values (antidictionary order)
varValues = fliplr(sort(varValues));
numOfValues = size(varValues, 2);

%check if index is valid
if 1<= index && index <= numOfValues^numOfSites
    
    pos = zeros(1, numOfSites);
    encoding = dec2base(index-1, numOfValues, numOfSites);
    for i=1:length(encoding)
        pos(i) = str2double(encoding(i));
    end
    
    pos = pos+1;
    
    %create the state
    state = '';
    for k=1:numOfSites 
        state = strcat(state, varValues{pos(k)});
    end

else
    error('the index is not in a valid stateSet range');
end

end