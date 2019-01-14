function D = myCellDivision(s)
%Cell Division Process Matrix
%HAVE TO BE CHECKED IF NOT SYMMETRIC CASE
%     UU   UM     UH     UF     MU     MM     MH     MF     HU     HM      HH   HF      FU      FM      FH      FF  
D =  [1    0      0      0      0      0      0      0      0      0       0    0       0       0       0       0;  %UU
      s    1-s    0      0      0      0      0      0      0      0       0    0       0       0       0       0;  %UM
      1-s  0      s      0      0      0      0      0      0      0       0    0       0       0       0       0;  %UH
      s    0      0      1-s    0      0      0      0      0      0       0    0       0       0       0       0;  %UF                                                                       
      s    0      0      0      1-s    0      0      0      0      0       0    0       0       0       0       0;  %MU
      0    1-s    0      0      s      0      0      0      0      0       0    0       0       0       0       0;  %MM
      0    0      s      0      1-s    0      0      0      0      0       0    0       0       0       0       0;  %MH
      0    0      0      s      1-s    0      0      0      0      0       0    0       0       0       0       0;  %MF                                    
      1-s  0      0      0      0      0      0      0      s      0       0    0       0       0       0       0;  %HU
      0    1-s    0      0      0      0      0      0      s      0       0    0       0       0       0       0;  %HM
      0    0      s      0      0      0      0      0      1-s    0       0    0       0       0       0       0;  %HH
      0    0      0      s      0      0      0      0      1-s    0       0    0       0       0       0       0;  %HF                                                   
      s    0      0      0      0      0      0      0      0      0       0    0       1-s     0       0       0;  %FU                                                                
      0    s      0      0      0      0      0      0      0      0       0    0       1-s     0       0       0;  %FM                                                                
      0    0      s      0      0      0      0      0      0      0       0    0       1-s     0       0       0;  %FH                                                                 
      0    0      0      s      0      0      0      0      0      0       0    0       1-s     0       0       0]; %FF                                                                



end

