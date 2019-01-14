function D = cellDivisionNew(s)
%Cell Division Process Matrix
%HAVE TO BE CHECKED IF NOT SYMMETRIC CASE
%     UU   UM     UH     MU     MM     MH     HU     HM     HH  
D =  [1    0      0      0      0      0      0      0      0        ;  %UU
      s    1-s    0      0      0      0      0      0      0        ;  %UM
      1-s  0      s      0      0      0      0      0      0        ;  %UH
      s    0      0      1-s    0      0      0      0      0        ;  %MU
      0    1-s    0      s      0      0      0      0      0        ;  %MM
      0    0      s      1-s    0      0      0      0      0        ;  %MH
      1-s  0      0      0      0      0      s      0      0        ;  %HU
      0    1-s    0      0      0      0      s      0      0        ;  %HM
      0    0      s      0      0      0      1-s    0      0       ;]; %HH



end

