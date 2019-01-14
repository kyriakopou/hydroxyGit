function T_in = createTotalObsPerDayTable(T_in)
%create extra columns in T table
%MAYBE THESE SHOULD GO TO A NEW TABLE NOT TO EDIT T
%BEC I GUESS IT HAS TO BE COPIED IN THIS CASE
T_in.bs0 = T_in.CC_bs0 + T_in.CT_bs0 + T_in.TC_bs0 + T_in.TT_bs0;
T_in.ox0 = T_in.CC_ox0 + T_in.CT_ox0 + T_in.TC_ox0 + T_in.TT_ox0;

T_in.bs3 = T_in.CC_bs3 + T_in.CT_bs3 + T_in.TC_bs3 + T_in.TT_bs3;
T_in.ox3 = T_in.CC_ox3 + T_in.CT_ox3 + T_in.TC_ox3 + T_in.TT_ox3;
% T_in.bs4 = T_in.CC_bs4 + T_in.CT_bs4 + T_in.TC_bs4 + T_in.TT_bs4;


T_in.bs6 = T_in.CC_bs6 + T_in.CT_bs6 + T_in.TC_bs6 + T_in.TT_bs6;
T_in.ox6 = T_in.CC_ox6 + T_in.CT_ox6 + T_in.TC_ox6 + T_in.TT_ox6;
% T_in.bs7 = T_in.CC_bs7 + T_in.CT_bs7 + T_in.TC_bs7 + T_in.TT_bs7;


end