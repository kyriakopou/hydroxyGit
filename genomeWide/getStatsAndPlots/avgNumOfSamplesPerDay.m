function avgObs = avgNumOfSamplesPerDay(T_in_ext)
%get the extended T_in table containing the total
%num of observations per day and compute the avg num of samples per day (bs, oxBS)

%total avg num of BS observations per day
avgObsBisPerDay1 = nanmean(T_in_ext.bs0);
avgObsBisPerDay2 = nanmean(T_in_ext.bs3);
avgObsBisPerDay3 = nanmean(T_in_ext.bs6);
%create a matrix for all days
avgObsBis = [avgObsBisPerDay1; avgObsBisPerDay2; avgObsBisPerDay3];
%total avg num of oxBS observations per day
avgObsOxPerDay1 = nanmean(T_in_ext.ox0);
avgObsOxPerDay2 = nanmean(T_in_ext.ox3);
avgObsOxPerDay3 = nanmean(T_in_ext.ox6);
%create matrix for all days 
avgObsOx = [avgObsOxPerDay1; avgObsOxPerDay2; avgObsOxPerDay3];

% avgObsOx = [];

avgObs = [avgObsBis, avgObsOx];


end