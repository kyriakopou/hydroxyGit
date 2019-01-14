function D = cellDivision(s)
%Cell Division Process Matrix

%     UU   UM     MU     UH     HU     HM     MH     MM     HH  
D =  [1    0      0      0      0      0      0      0      0        ;  %UU
      s    1-s    0      0      0      0      0      0      0        ;  %UM
      1-s  0      s      0      0      0      0      0      0        ;  %MU
      s    0      0      1-s    0      0      0      0      0        ;  %UH
      1-s  0      0      0      s      0      0      0      0        ;  %HU
      0    1-s    0      0      s      0      0      0      0        ;  %HM
      0    0      s      1-s    0      0      0      0      0        ;  %MH
      0    1-s    s      0      0      0      0      0      0        ;  %MM
      0    0      0      1-s    s      0      0      0      0       ;]; %HH



end

