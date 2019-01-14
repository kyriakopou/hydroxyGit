function numResultsRow = getNumerResults(params, pAllStates, maxObsDaysVector)

        %store hidden state probabilities properly for the output
        %pAllStates has here always 3 rows (days without observ have null entries)
        uuCol = pAllStates(:,1)';
        umCol = sum(pAllStates(:,2:3), 2)';
        tothCol = sum(pAllStates(:,4:5), 2)' + sum(pAllStates(:,6:7), 2)' + pAllStates(:,9)';
        mmCol = pAllStates(:,8)';

        %store the parameters (efficiencies) values
        maint_b0 = params(1);
        maint_b1 = params(2);
        deNovo_b0 = params(3);
        deNovo_b1 = params(4);
        hydroxy_b0 = params(5);
        hydroxy_b1 = params(6);
        pRecogn = params(7);
               
        TimePointsArray = repmat(maxObsDaysVector', 1);
        maintEffic = repmat(maint_b0, 1, 3) + repmat(maint_b1, 1, 3) .* TimePointsArray;
        deNovoEffic = repmat(deNovo_b0, 1, 3) + repmat(deNovo_b1, 1, 3) .* TimePointsArray;
        hydroxyEffic = repmat(hydroxy_b0, 1, 3) + repmat(hydroxy_b1, 1, 3) .* TimePointsArray;
        
        %store results into table's row
        numResultsRow = [mmCol(:,1), tothCol(:,1), umCol(:,1), uuCol(:,1), ...
            mmCol(:,2), tothCol(:,2), umCol(:,2), uuCol(:,2), mmCol(:,3), tothCol(:,3), umCol(:,3), uuCol(:,3), ...
            maintEffic(:,1), deNovoEffic(:,1), hydroxyEffic(:,1), maintEffic(:,2), deNovoEffic(:,2), hydroxyEffic(:,2), ...
            maintEffic(:,3), deNovoEffic(:,3), hydroxyEffic(:,3), pRecogn];
        

end