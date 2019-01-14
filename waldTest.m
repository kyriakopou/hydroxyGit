function [paramNames, ] = waldTest(params, testCov)

%define the g function and the Jacobian for the rate params of interest (not-nuisance)
    Wald = zeros(1, floor(size(params, 2) / 2));
    for testParam = 2:2:size(params, 2)
        g = params(testParam);
        Jacob = zeros(1, size(params,2));
        Jacob(testParam) = 1;
        Wald(testParam / 2) = g * inv(Jacob * CovOut * Jacob')*g';
        %convert cell array position to string (char array)
        paramNames{testParam / 2} = char(paramNames(testParam));
        fprintf(fileID, 'Wald.Stat(%s) = %d \n', paramNames, Wald(testParam / 2));
    end

    %wald test for non-recognition prob. if it is in the keeped params
    if (ub(end) ~= lb(end)) 
        gP = params(end) - 1;
        JacobP = zeros(1, size(params, 2));
        JacobP(7) = 1;
        WaldP = gP * inv(JacobP*CovOut*JacobP')*gP';
        fprintf(fileID, 'Wald.Stat(p) = %d\n', WaldP);
    end

    %wald test for total methylation lambda
    gLambda = [params(2) + params(4) - params(2)*params(3) - params(1)*params(4);
                params(2)*params(4)];
    JacobLambda = [-params(4), 1-params(3), -params(2), 1-params(1), zeros(1, size(params, 2)-4);
                   0, params(4), 0, params(2), zeros(1, size(params, 2) - 4)];
    WaldLambda = gLambda' * inv(JacobLambda*Cov*JacobLambda')*gLambda;


end