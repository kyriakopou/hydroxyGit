function [T, notNanRows] = filteringTable(T, minDayCoverage)

T_final_ext = createTotalObsPerDayTable(T);
    
day0ind = (T_final_ext.bs0 < minDayCoverage | T_final_ext.ox0 < minDayCoverage);
day3ind = (T_final_ext.bs3 < minDayCoverage | T_final_ext.ox3 < minDayCoverage);
day6ind = (T_final_ext.bs6 < minDayCoverage | T_final_ext.ox6 < minDayCoverage);

%make nan day0 BS
T.CC_bs0(day0ind) = nan; 
T.CT_bs0(day0ind) = nan;
T.TC_bs0(day0ind) = nan;
T.TT_bs0(day0ind) = nan;
%make nan day0 oxBS
T.CC_ox0(day0ind) = nan; 
T.CT_ox0(day0ind) = nan;
T.TC_ox0(day0ind) = nan;
T.TT_ox0(day0ind) = nan;
%make nan day3 BS
T.CC_bs3(day3ind) = nan; 
T.CT_bs3(day3ind) = nan;
T.TC_bs3(day3ind) = nan;
T.TT_bs3(day3ind) = nan;
%make nan day3 oxBS
T.CC_ox3(day3ind) = nan; 
T.CT_ox3(day3ind) = nan;
T.TC_ox3(day3ind) = nan;
T.TT_ox3(day3ind) = nan;
%make nan day6 BS
T.CC_bs6(day6ind) = nan; 
T.CT_bs6(day6ind) = nan;
T.TC_bs6(day6ind) = nan;
T.TT_bs6(day6ind) = nan;
%make nan day6 oxBS
T.CC_ox6(day6ind) = nan; 
T.CT_ox6(day6ind) = nan;
T.TC_ox6(day6ind) = nan;
T.TT_ox6(day6ind) = nan;

%delete all nan rows
notNanRows = all(~isnan(T{:, 3:end}), 2);
T = T(notNanRows, :);


end