%we create a cell array containing the states of the stateSet as strings
%(WORKS FOR TWO SITES - HAS TO BE EXTENDED FOR MORE SITES)
function stateSet = createStateSet(varValues, numOfSites)

numOfValues = size(varValues, 2);
stateSet = cell(1, numOfValues^numOfSites);

for i=1:numOfValues
    for j=1:numOfValues
        stateSet{(i-1)*numOfValues + j} = strcat(varValues{i}, varValues{j});
    end   
end


end