function processAll = getProcessForEachDay(days, process)
%modify process vector to get entries 
%for all dataPoints

processAll = cell(days(end), 1);
numOfDays = length(days);

for i=2:numOfDays
    processAll{days(i)} = process(i);
end

for i=days(end)-1:-1:1
    if (isempty(processAll{i}))
        processAll{i} = processAll{i+1};
    end
end    



end