function [oneDayRows, twoDaysRows, threeDaysRows] = getRowsOfTable(T)

%distinguish between input and output file cases
if strcmp(T.Properties.VariableNames{2}, 'position') 
    firstDayRows = ~isnan(T.CC_bs0);
    secDayRows = ~isnan(T.CC_bs3);
    thirdDayRows = ~isnan(T.CC_bs6);
else
    firstDayRows = ~isnan(T.mm_d0);
    secDayRows = ~isnan(T.mm_d3);
    thirdDayRows = ~isnan(T.mm_d6);
end


oneDayRows = firstDayRows + secDayRows + thirdDayRows == 1;
twoDaysRows = firstDayRows + secDayRows + thirdDayRows == 2;
threeDaysRows = firstDayRows + secDayRows + thirdDayRows == 3;

end