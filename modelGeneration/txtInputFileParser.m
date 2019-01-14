function [model, vars, params, maxDegree, knots, process, rulesLHS, rulesRHS, probs] = txtInputFileParser(fileName)

fid = fopen(fileName);
%read the number of processes
tline = fgetl(fid);
numOfProc = str2double(tline(6));
%initialize output matrices
process = cell(1, numOfProc);
rulesLHS = cell(1, numOfProc);
rulesRHS = cell(1, numOfProc);
probs = cell(1, numOfProc);
%read the type of the model
tline = fgetl(fid);
a = strsplit(tline, '=');
model = a{2};
%read model's vars for the sites
v = strsplit(fgetl(fid), {'=', ','});
vars = v(2:end);
%read model's params names
p = strsplit(fgetl(fid), {'=', ','});
params = p(2:end);
%read model's maxDegree
m = strsplit(fgetl(fid), {'='});
maxDegree = str2double(m(2));
%read the knots
k = strsplit(fgetl(fid), {'=', '[', ']' ','});
knots = str2double(k(2:end-1));

i=1;
fgetl(fid);
tline = fgetl(fid);
while i <= numOfProc
   process{i} = tline; 
   %read first rule 
   tline = fgetl(fid);
   while ~strcmp(tline, '---') && ~strcmp(tline, '--end--')
       %split the string according to specified delimeters
       a = strsplit(tline, {'->', ' @'});
       
       rulesLHS{i} = [rulesLHS{i}, a(1)];
       rulesRHS{i} = [rulesRHS{i}, a(2)];
       probs{i} = [probs{i}, a(3)];
       tline = fgetl(fid);
   end
   i = i+1;
   tline = fgetl(fid);
end


%make now the checks to the file
if strcmp(model, 'dtmc')
    for i=1:numOfProc
        uniqueRulesLHS = unique(rulesLHS{i});
        for u=uniqueRulesLHS
            if sum(sym(probs{i}(strcmp(rulesLHS{i}, u)))) ~= 1
                error('the sum of transition probabilities from a is not one')
            end
        end
    end
end
fclose(fid);

end