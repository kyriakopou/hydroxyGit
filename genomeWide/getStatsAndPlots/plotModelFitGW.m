function plotModelFitGW(T_in, T_out, EMatrix, daysLabels)

%FOR THE WHOLE FIT SECTION WE HAVE TO CHANGE THE MODEL RESULTS BY
%CONSIDERING THE ERRORS (WE NEED FOR THIS ALL HIDDENS STATES - NEW RUN)
%--FIT DATA VS MODEL PLOTS-- 


numOfDays = length(daysLabels);
%last index is the day
pAllStates = zeros(size(T_out, 1), 9, numOfDays);
pAllStates(:,:,1) = T_out{:,5:13};
pAllStates(:,:,2) = T_out{:,14:22};
pAllStates(:,:,3) = T_out{:,23:31};

%compute pModel
pModel = zeros(size(T_out, 1), 4, numOfDays);
for k=1:numOfDays
    for e=1:2
        pModel(:,:,k,e) = pAllStates(:,:,k) * EMatrix(:,:,k,e);
    end
end


%create extra columns in T table
%MAYBE THESE SHOULD GO TO A NEW TABLE NOT TO EDIT T
%BEC I GUESS IT HAS TO BE COPIED IN THIS CASE
T_in.bs0 = T_in.CC_bs0 + T_in.CT_bs0 + T_in.TC_bs0 + T_in.TT_bs0;
T_in.ox0 = T_in.CC_ox0 + T_in.CT_ox0 + T_in.TC_ox0 + T_in.TT_ox0;

T_in.bs3 = T_in.CC_bs3 + T_in.CT_bs3 + T_in.TC_bs3 + T_in.TT_bs3;
T_in.ox3 = T_in.CC_ox3 + T_in.CT_ox3 + T_in.TC_ox3 + T_in.TT_ox3;

T_in.bs6 = T_in.CC_bs6 + T_in.CT_bs6 + T_in.TC_bs6 + T_in.TT_bs6;
T_in.ox6 = T_in.CC_ox6 + T_in.CT_ox6 + T_in.TC_ox6 + T_in.TT_ox6;


%%
%for the hemimethyl bar plots we can not compare data and output since we
%need MLE or normalization of frequencies of both experiments
TT_bsd0 = T_in.TT_bs0 ./ T_in.bs0;
TT_bsd3 = T_in.TT_bs3 ./  T_in.bs3;
TT_bsd6 = T_in.TT_bs6 ./ T_in.bs6;

TT_oxd0 = T_in.TT_ox0 ./ T_in.ox0;
TT_oxd3 = T_in.TT_ox3 ./  T_in.ox3;
TT_oxd6 = T_in.TT_ox6 ./ T_in.ox6;

figure();
subplot(3,2,1);
TT_bsAll = {[TT_bsd0, TT_bsd3, TT_bsd6]; ...
    pModel(:,1,1:3,1)};
aboxplot(TT_bsAll, 'colorgrad', 'blue_up', 'labels', daysLabels);
legend('data', 'model'); % Add a legend
title('TT (BS)')

subplot(3,2,2);
TT_oxAll = {[TT_oxd0, TT_oxd3, TT_oxd6]; ...
    pModel(:,1,1:3,2)};
aboxplot(TT_oxAll, 'colorgrad', 'blue_up', 'labels', daysLabels);
legend('data', 'model'); % Add a legend
title('TT (oxBS)')

%%
%for the hemimethyl bar plots we can not compare data and output since we
%need MLE or normalization of frequencies of both experiments
hemi_bsd0 = (T_in.CT_bs0 + T_in.TC_bs0) ./ T_in.bs0;
hemi_bsd3 = (T_in.CT_bs3 + T_in.TC_bs3) ./  T_in.bs3;
hemi_bsd6 = (T_in.CT_bs6 + T_in.TC_bs6) ./ T_in.bs6;

hemi_oxd0 = (T_in.CT_ox0 + T_in.TC_ox0) ./ T_in.ox0;
hemi_oxd3 = (T_in.CT_ox3 + T_in.TC_ox3) ./  T_in.ox3;
hemi_oxd6 = (T_in.CT_ox6 + T_in.TC_ox6) ./ T_in.ox6;

% figure();
subplot(3,2,3);
hemi_bsAll = {[hemi_bsd0, hemi_bsd3, hemi_bsd6]; ...
    pModel(:,2,1:3,1) + pModel(:,3,1:3,1)};
aboxplot(hemi_bsAll, 'colorgrad', 'green_up', 'labels', daysLabels);
legend('data', 'model'); % Add a legend
title('TC+CT (BS)')


% figure();
subplot(3,2,4);
hemi_oxAll = {[hemi_oxd0, hemi_oxd3, hemi_oxd6]; ...
    pModel(:,2,1:3,2) + pModel(:,3,1:3,2)};
aboxplot(hemi_oxAll, 'colorgrad', 'green_up', 'labels', daysLabels);
legend('data', 'model'); % Add a legend
title('TC+CT (oxBS)')


%%
%CC bs
CC_bs0 = T_in.CC_bs0 ./ T_in.bs0;
CC_bs3 = T_in.CC_bs3 ./ T_in.bs3;
CC_bs6 = T_in.CC_bs6 ./ T_in.bs6;
%CC ox
CC_ox0 = T_in.CC_ox0 ./ T_in.ox0;
CC_ox3 = T_in.CC_ox3 ./ T_in.ox3;
CC_ox6 = T_in.CC_ox6 ./ T_in.ox6;


subplot(3,2,5);
CC_bsAll = {[CC_bs0, CC_bs3, CC_bs6]; ...
    pModel(:,4,1:3,1)};
aboxplot(CC_bsAll, 'colorgrad', 'red_up', 'labels', daysLabels);
legend('data', 'model'); % Add a legend
title('CC (BS)')

subplot(3,2,6);
CC_oxAll = {[CC_ox0, CC_ox3, CC_ox6]; ...
    pModel(:,4,1:3,2)};
aboxplot(CC_oxAll, 'colorgrad', 'red_up', 'labels', daysLabels);
legend('data', 'model'); % Add a legend
title('CC (oxBS)')




end

