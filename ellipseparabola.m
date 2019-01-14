function [c, ceq, gradc, gradceq] = ellipseparabola(x, maxDegree, knots, ubR, t_max)

numOfParams = length(x);
numOfKnots = size(knots, 2);
numOfCoef = maxDegree + 1 + numOfKnots;
numOfRates = (numOfParams-1) / numOfCoef;


%define efficiency f, df_dt, and df_dt_dt as symbolic expression and find tArgMinMax
syms t
b = sym('b', [1 numOfCoef]);

polyTerms = t.^(0:maxDegree);
truncTerms = sym('trunc', [1, numOfKnots]);
for i=1:numOfKnots
    truncTerms(i) = piecewise((t >= knots(i)), b(maxDegree+1+i)*(t - knots(i)).^(maxDegree), t < knots(i), 0);
end

sum(truncTerms);
f = sum(b(1:maxDegree+1) .* polyTerms) + sum(truncTerms); 

% f = sum(b .* baseTerms);

df_dt = diff(f, t);
df_dt_dt = diff(df_dt, t);


% numOfConstraints = numOfRates * 2*(numOfKnots+1);
numOfConstraints = 2*numOfRates;

c = zeros(numOfConstraints, 1);
ceq = [];
gradc = zeros(numOfParams, numOfConstraints);
gradceq = [];

for i=1:numOfRates
    
    constBase = (i-1)*2;
    coefBase = (i-1)*numOfCoef;
    
    rateParams = x((i-1)*numOfCoef+1:i*numOfCoef);
    
    %has to be generalized for arbitrary numOfKnots
    if ~isempty(knots)
        knotsMatrix = [knots(2) t_max; knots(1) knots(2); 0 knots(1)];
    else
        knotsMatrix = [0 t_max];
    end
        
    %for the function to be printed
    f_b = subs(f, b, rateParams);
    
%     figure(1);
%     ezplot(f_b, [0 12 -1 4]);
    
    %symbolic solutions of the derivative (general form with 3 different possible solutions 
    %- some of them might be invalid)
    %evaluating first f at beta points would not work because of numerical
    %issues
%     [solx, param, cond] = solve(df_dt == 0, t, 'ReturnConditions', true);
        
    tArgMinMax = solve(df_dt == 0, t);

    %assign only to terms with non zero denominator
    [num, dem] = numden(tArgMinMax);
    demEval = double(subs(dem, b, rateParams));
    ind = (demEval ~= 0);
    
    %consider only the points where der = 0 solution is valid
    tArgMinMax = tArgMinMax(ind);
    knotsMatrix = knotsMatrix(ind, :);
    tArgMinMax_b = double(subs(tArgMinMax, b, rateParams));
    %consider only the points which are within the right knotBounds
    if ~isempty(tArgMinMax)
        tArgMinMax_b = double(subs(tArgMinMax, b, rateParams));
        indBounds = (knotsMatrix(:,1) <= tArgMinMax_b & tArgMinMax_b <= knotsMatrix(:,2));
        tArgMinMax = tArgMinMax(indBounds);
        tArgMinMax_b = tArgMinMax_b(indBounds);
    end
    
    %add the first and the last point
    criticalPoints = [tArgMinMax_b; 0; t_max];
    %and to the output of f for 
    f_criticalPoints = subs(f, t, criticalPoints);
    f_criticalPoints_b = subs(f_b, t, criticalPoints);
                     
    %compute symb expression of df_dt_dt at t_ArgMinMax 
    df_dt_dt_tArgMinMax = subs(df_dt_dt, t, tArgMinMax_b);
    df_dt_dt_tArgMinMax_b = subs(df_dt_dt_tArgMinMax, b, rateParams);

    %remove nans (points where sec der is not defined or it is zero)
%         indNotNanNotZero = (~isnan(df_dt_dt_tArgMinMax_b) & df_dt_dt_tArgMinMax_b ~= 0);
%         df_dt_dt_tArgMinMax_b = df_dt_dt_tArgMinMax_b(indNotNanNotZero);
% 
%         f_tArgMinMax = f_tArgMinMax(indNotNanNotZero);
%         tArgMinMax_b = tArgMinMax_b(indNotNanNotZero);

    %evaluate sec derivative at non-nan non-zero tArgMinMax points

%     figure()
%     ezplot(fEval_b, [0 12 -1 4]);

    %if there is tArgMinMax entry within the bounds get the
    %conditions and the derivatives        
    %do we have min max or both?
    [fMax, indMax] = max(f_criticalPoints_b);
    [fMin, indMin] = min(f_criticalPoints_b);

%             if df_dt_dt_tArgMinMax_b(e) < 0
    c(constBase+1) = fMax - ubR;
%             else 
    c(constBase+2) = -fMin;
%             end

    %derivatives of the nonlinear constraints
    gradsMax = zeros(numOfCoef, 1);
    gradsMin = zeros(numOfCoef, 1);
    for k=1:numOfCoef 
        df_critical_max_db = diff(f_criticalPoints(indMax), b(k));
        gradsMax(k) = subs(df_critical_max_db, b, rateParams); 
        %this has to be checked 
        df_critical_min_db = diff(f_criticalPoints(indMin), b(k));
        gradsMin(k) = subs(df_critical_min_db, b, rateParams);
    end
    %max
%       if df_dt_dt_tArgMinMax_b(e) < 0
    gradc(:, constBase+1) = [zeros(coefBase,1); gradsMax; zeros(numOfParams-(numOfCoef+coefBase), 1)]; 
    %min    
%       else                   
    gradc(:, constBase+2) = [zeros(coefBase,1); -gradsMin; zeros(numOfParams-(numOfCoef+coefBase), 1)];
%       end
        
        
end


end