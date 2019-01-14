function [M1, dermu2M1, dermu1M1, derhM1] = maintenance(mu1)
%Maintenace Methylation Process Matrix
%mu1 is maintenance methylation probability


%     UU        UM             MU          UH       HU       HM       MH        MM        HH  
M1 = [1         0              0           0        0        0        0         0         0;  %UU
      0         1-mu1          0           0        0        0        0         mu1       0;  %UM
      0         0              1-mu1       0        0        0        0         mu1       0;  %MU
      0         0              0           1        0        0        0         0         0;  %UH
      0         0              0           0        1        0        0         0         0;  %HU
      0         0              0           0        0        1        0         0         0;  %HM
      0         0              0           0        0        0        1         0         0;  %MH
      0         0              0           0        0        0        0         1         0;  %MM
      0         0              0           0        0        0        0         0         1;];%HH 
  
  
%          UU         UM             MU          UH       HU       HM       MH       MM        HH         
dermu1M1 = [0         0              0           0        0        0        0        0         0;  %UU
            0         -1             0           0        0        0        0        1         0;  %UM
            0         0              -1          0        0        0        0        1         0;  %MU
            0         0              0           0        0        0        0        0         0;  %UH
            0         0              0           0        0        0        0        0         0;  %HU
            0         0              0           0        0        0        0        0         0;  %HM
            0         0              0           0        0        0        0        0         0;  %MH
            0         0              0           0        0        0        0        0         0;  %MM
            0         0              0           0        0        0        0        0         0;];%HH     
               
      
dermu2M1 = zeros(9,9);              
derhM1 = zeros(9,9);   
        

end