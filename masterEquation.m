function dp_dv = masterEquation(time, p_v, x, maxDegree, knots, numOfHiddenStates, derCompFlag, QMatrixFunction)

numOfParams = length(x);
%compute the efficiencies at timepoint time
EFF = getEfficiencies(x, maxDegree, time, knots);

mu2 = EFF(:,2);
h = EFF(:,3);
f = EFF(:,4);
dem = EFF(:,5);

%get Q from the time varying rates
[Q, derQ, secDerQ] = QMatrixFunction(dem, f, h, mu2, time);

%split the input into the transient prob vector (first 16 entries)
p = p_v(1:numOfHiddenStates);
dp_dv = Q'*p;

firstDer = p_v(numOfHiddenStates+1:(numOfParams+1)*numOfHiddenStates);

if derCompFlag >= 1
     for i=1:numOfParams
        %derivative vector for each param
        derpI = firstDer(numOfHiddenStates*(i-1)+1:i*numOfHiddenStates);
        %get new prob vector and derivative vector for each parameter
        dp_dv = [dp_dv; Q'*derpI + derQ(:,:,i)'*p];
          
    end
    if derCompFlag == 2
        secDer = p_v((numOfParams+1)*numOfHiddenStates+1:end);
        for i=1:numOfParams
            for j=1:numOfParams
                %second derivative vector
                secDerpIJ = secDer((i-1)*numOfParams+(j-1)*numOfHiddenStates+1:(i-1)*numOfParams+j*numOfHiddenStates);
                derpI = firstDer(numOfHiddenStates*(i-1)+1:i*numOfHiddenStates);
                derpJ = firstDer(numOfHiddenStates*(j-1)+1:j*numOfHiddenStates);

                dp_dv = [dp_dv; Q'* secDerpIJ + derQ(:,:,i)' * derpJ + derQ(:,:,j)' * derpI + secDerQ(:,:,i,j) * p];
            end    
        end
    end    
end

end