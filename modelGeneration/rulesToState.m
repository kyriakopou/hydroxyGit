function [descendants, probSym] = rulesToState(state, varLhs, varRhs, prob)
%compute the rules that are possible from state
%s and compute the descendants and the rates or the prob to each of them.


%make the controls
if size(varLhs, 2) ~= size(varRhs, 2) || size(varLhs, 2) ~= size(prob, 2)
    error('varLhs, varRhs and prob cell arrays are not of the same size');
end

%find the entry of cell array varLhs that can be applied to state
%(WE ASSUME ONLY ONE RULE APPLIES TO EACH STATE--NEEDS TO BE GENERALIZED)
numOfRules = size(varLhs, 2);
indMatches = cell(1, numOfRules);
%for a rule that contains at least a dyad in lhs
rulesChangeVarInd = cell(1, numOfRules);

for i=1:numOfRules
    %if we find a substring of varLhs{i} inside state's string keep rule i
    %indMatches is the the first index of each occurence of rule i in state
    indMatches{i} = strfind(state, varLhs{i});
    %rulesChangesVarInd is the index of the element of match that can be
    %changed (e.g. 1st or 2nd element for um, mu case)
    %THE REACTION THAT CHANGES THE VARIABLE HAS TO BE FIRST FOR THE
    %FOLLOWING TO WORK
    if ~isempty(find(varLhs{i}~=varRhs{i}))
        rulesChangeVarInd{i} = find(varLhs{i}~=varRhs{i});
    else
        rulesChangeVarInd{i} = rulesChangeVarInd{i-1};
    end    
end    

%keep only all reactions that can be applied in the state
nonEmptyInd = find(~cellfun(@isempty, indMatches));
indMatches = indMatches(:, nonEmptyInd);
rulesChangeVarInd = rulesChangeVarInd(:, nonEmptyInd);
varLhs = varLhs(:, nonEmptyInd);
varRhs = varRhs(:, nonEmptyInd);
prob = prob(:, nonEmptyInd);


%if there is at least one matching rule
if ~isempty(nonEmptyInd)

        %compute with how many ways one lhs can be changed
        %find varLhS unique values
        unVarLhs = fliplr(unique(varLhs));
        %number of matches for each unique lhs
        numOfMatches = zeros(1, size(unVarLhs, 2));
        %numberOfRulesLhs{k} equals also the number of possible
        %modifications for lhs{k}
        numOfRulesForLhs = cell(1, size(unVarLhs, 2)); 
        statesChangedVarsIndices = cell(1, size(unVarLhs, 2));
        numOfDescendants = zeros(1, size(unVarLhs, 2));
        for k=1:size(unVarLhs, 2)
            %indices of matches for lhs k
            statesChangedVarsIndices{k} = indMatches{:, find(strcmp(varLhs, unVarLhs{k}), 1)} + rulesChangeVarInd{find(strcmp(varLhs, unVarLhs{k}), 1)} -1;
            
            %numOf possible transformations for lhs k
            numOfRulesForLhs{k} = size(varLhs(strcmp(varLhs, unVarLhs{k})), 2);
            numOfMatches(k) = size(indMatches{:, find(strcmp(varLhs, unVarLhs{k}), 1)}, 2);
            
            %numOf possible Descendants of the state due to change of lhs k
            numOfDescendants(k) = numOfRulesForLhs{k}^numOfMatches(k);
                   
        end
        
        newStatesChangedVars = cell(1, size(unVarLhs, 2));
%         statesChangedVarsIndices = cell(1, size(unVarLhs, 2));
        totNumOfDescendants = prod(numOfDescendants);
        descendants = cell(1, totNumOfDescendants);
        probSym = sym(zeros(1, totNumOfDescendants));
        
        %----HARDCODED----
        %WORKS ONLY FOR ONE or TWO LHS VARIABLES--NEEDS TO BE GENERALIZED (FOR MY RULES NOW THIS IS ALWAYS ENOUGH) 
        if size(unVarLhs, 2) == 2
            for i=1:numOfDescendants(1)
                rhsChangeVars = varRhs(strcmp(varLhs, unVarLhs{1}));
                rhsChangeVars = cellfun(@(x) x(rulesChangeVarInd{find(strcmp(varLhs, unVarLhs{1}), 1)}), rhsChangeVars);
                newStatesChangedVars{1} = getStateFromIndex(i, cellstr(rhsChangeVars')', numOfMatches(1));
                           
                %the probability for this transition
                newStatesChangedVars{1} = reshape(newStatesChangedVars{1}, [], 1);
                transProd1 = '';
                for k=1:length(newStatesChangedVars{1})
                    %find the prob of this transition for newStatesUnfixedVars(k) to happen 
%                   pTemp = prob{strcmp(varLhs, unVarLhs{1}) & strcmp(varRhs, newStatesChangedVars{1}(k))};
                    changeVars = cellfun(@(x) x(rulesChangeVarInd{find(strcmp(varLhs, unVarLhs{1}), 1)}), varRhs(strcmp(varLhs, unVarLhs{1})));
                    ind = strfind(changeVars, newStatesChangedVars{1}(k));
                    pTemp = prob{ind}; 
                    
                    if isempty(transProd1) 
                        transProd1 = strcat('(', pTemp, ')');
                    else    
                        transProd1 = strcat(transProd1, '*', '(', pTemp, ')');
                    end    
                end
                
                for j=1:numOfDescendants(2)                      
                    %change (or not) the variables that the rule should be applied
                    rhsChangeVars = varRhs(strcmp(varLhs, unVarLhs{2}));
                    rhsChangeVars = cellfun(@(x) x(rulesChangeVarInd{find(strcmp(varLhs, unVarLhs{2}), 1)}), rhsChangeVars);
                    newStatesChangedVars{2} = getStateFromIndex(j, cellstr(rhsChangeVars')', numOfMatches(2));        

                    %compute the descendant
                    ind = (i-1)*(numOfDescendants(2)) + j;
                    descendants{ind} = state;
                    descendants{ind}(statesChangedVarsIndices{1}) = newStatesChangedVars{1};
                    descendants{ind}(statesChangedVarsIndices{2}) = newStatesChangedVars{2};
                    
                    
                    %the probability for this transition
                    newStatesChangedVars = reshape(newStatesChangedVars, [], 1);
                    transProd2 = '';
                    for k=1:length(newStatesChangedVars{2})
                        %find the prob of this transition for newStatesUnfixedVars(k) to happen 
                        %WORKAROUND FOR THE 2 RULES CASE AGAIN
                        %changeVars of the second unique rule
                        changeVars = cellfun(@(x) x(rulesChangeVarInd{find(strcmp(varLhs, unVarLhs{2}), 1)}), varRhs(strcmp(varLhs, unVarLhs{2})));
                        %index of the changedVariable that was drawn in
                        %newStatesChangedVars
                        indNSC = numOfDescendants(1) + strfind(changeVars, newStatesChangedVars{2}(k));
                        pTemp = prob{indNSC}; 
                        if isempty(transProd2) 
                            transProd2 = strcat('(', pTemp, ')');
                        else    
                            transProd2 = strcat(transProd2, '*', '(', pTemp, ')');
                        end    
                    end
                    
                    transProd = strcat(transProd1,'*', transProd2); 
                    probSym(ind) = evalin(symengine, transProd);
                end
            end
        else
            for i=1:numOfDescendants(1)          
                rhsChangeVars = varRhs(strcmp(varLhs, unVarLhs{1}));
                rhsChangeVars = cellfun(@(x) x(rulesChangeVarInd{1}), rhsChangeVars);
                newStatesChangedVars = getStateFromIndex(i, cellstr(rhsChangeVars')', numOfMatches(1));
                
                %compute the descendant
                descendants{i} = state;
                descendants{i}(statesChangedVarsIndices{1}) = newStatesChangedVars;

                %the probability for this transition
                newStatesChangedVars = reshape(newStatesChangedVars, [], 1);
                transProd = '';
                %for each of the matches
                for k=1:length(newStatesChangedVars)
                    %find the prob of this transition for newStatesUnfixedVars(k) to happen 
                    %FOLLOWING LINE IS ALSO WORKAROUND
                    changeVars = cellfun(@(x) x(rulesChangeVarInd{1}), varRhs(strcmp(varLhs, unVarLhs{1})));
                    ind = strfind(changeVars, newStatesChangedVars(k));
                    pTemp = prob{ind}; 
                    if isempty(transProd) 
                        transProd = strcat('(', pTemp, ')');
                    else    
                        transProd = strcat(transProd, '*', '(', pTemp, ')');
                    end    
                end

                probSym(i) = evalin(symengine, transProd);
                 
            end    
        end
        
    
else

    descendants = cell(0);
    probSym = sym('probs', [0 0]);

end



end