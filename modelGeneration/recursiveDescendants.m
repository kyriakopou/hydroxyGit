function [newStatesUnfixedVars] = recursiveDescendants(unLHS, numOfDescendants, rhsNew, numOfMatches)

% rhsNew, indMatches, numOfMatches, rulesChangeVarInd

if unLHS > 1
    for i=1:numOfDescendants
        newStatesUnfixedVars = getStateFromIndex(i, rhsNew{unLHS}, numOfMatches{unLHS});
        newStatesUnfixedVarsNew = recursiveDescendants(unLHS-1, numOfDescendants(unLHS-1), rhsNew{unLHS-1}, numOfMatches{unLHS-1});
        newStatesUnfixedVars = {newStatesUnfixedVars; newStatesUnfixedVarsNew};
        
    end                 
else
    for i=1:numOfDescendants
        newStatesUnfixedVars = getStateFromIndex(i, rhsNew{unLHS}, numOfMatches{unLHS});
    end    
        
end