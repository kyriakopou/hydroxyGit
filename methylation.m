function [M, derM, secondDerM] = methylation(mu1, mu2, p, t, derCompFlag)
%Methylation Process Matrix
%mu1 is maintenance methylation
%mu2 is de-Novo methylation
%p is the probability that h is recognized as u


tot = 2*mu2-mu2^2;   

pHemi = mu1+mu2-mu1*mu2;


%     UU        UM                  MU                UH                            HU                           HM                      MH                     MM                HH  
M  = [1-tot     mu2-mu2^2           mu2-mu2^2         0                             0                            0                       0                      mu2^2             0;  %UU
      0         1-pHemi             0                 0                             0                            0                       0                      pHemi             0;  %UM
      0         0                   1-pHemi           0                             0                            0                       0                      pHemi             0;  %MU
      0         0                   0                 p*(1-mu2)+(1-p)*(1-pHemi)     0                            0                       p*mu2+(1-p)*pHemi      0                 0;  %UH
      0         0                   0                 0                             p*(1-mu2)+(1-p)*(1-pHemi)    p*mu2+(1-p)*pHemi       0                      0                 0;  %HU
      0         0                   0                 0                             0                            1                       0                      0                 0;  %HM
      0         0                   0                 0                             0                            0                       1                      0                 0;  %MH
      0         0                   0                 0                             0                            0                       0                      1                 0;  %MM
      0         0                   0                 0                             0                            0                       0                      0                 1;];%HH 
 
derM = [];
secondDerM = [];

if derCompFlag  
    %            UU         UM             MU          UH             HU             HM             MH            MM        HH         
    derb0mu1M  = [0         0              0           0              0              0              0             0         0;  %UU
                  0         -1+mu2         0           0              0              0              0             1-mu2     0;  %UM
                  0         0              -1+mu2      0              0              0              0             1-mu2     0;  %MU
                  0         0              0           (-1+p)*(1-mu2) 0              0              (1-p)*(1-mu2) 0         0;  %UH
                  0         0              0           0              (-1+p)*(1-mu2) (1-p)*(1-mu2)  0             0         0;  %HU
                  0         0              0           0              0              0              0             0         0;  %HM
                  0         0              0           0              0              0              0             0         0;  %MH
                  0         0              0           0              0              0              0             0         0;  %MM
                  0         0              0           0              0              0              0             0         0;];%HH     

    %             UU        UM             MU          UH                HU            HM          MH            MM        HH         
    derb0mu2M =  [2*mu2-2   1-2*mu2        1-2*mu2     0                 0             0           0             2*mu2     0;  %UU
                  0         -1+mu1         0           0                 0             0           0             1-mu1     0;  %UM
                  0         0              -1+mu1      0                 0             0           0             1-mu1     0;  %MU
                  0         0              0           -1+mu1*(1-p)      0             0           1-mu1*(1-p)   0         0;  %UH
                  0         0              0           0                 -1+mu1*(1-p)  1-mu1*(1-p) 0             0         0;  %HU
                  0         0              0           0                 0             0           0             0         0;  %HM
                  0         0              0           0                 0             0           0             0         0;  %MH
                  0         0              0           0                 0             0           0             0         0;  %MM
                  0         0              0           0                 0             0           0             0         0;];%HH   


    %         UU        UM             MU          UH           HU          HM           MH            MM        HH         
    derpM  = [0         0              0           0            0           0            0             0         0;  %UU
              0         0              0           0            0           0            0             0         0;  %UM
              0         0              0           0            0           0            0             0         0;  %MU
              0         0              0           mu1*(1-mu2)  0           0            -mu1*(1-mu2)  0         0;  %UH
              0         0              0           0            mu1*(1-mu2) -mu1*(1-mu2) 0             0         0;  %HU
              0         0              0           0            0           0            0             0         0;  %HM
              0         0              0           0            0           0            0             0         0;  %MH
              0         0              0           0            0           0            0             0         0;  %MM
              0         0              0           0            0           0            0             0         0;];%HH 


    derb1mu1M = derb0mu1M * t;
    derb1mu2M = derb0mu2M * t;    
    derb0hM = zeros(9, 9);
    derb1hM = zeros(9, 9); 


    %concatenate  all derivative matrices into a 3 dimensional matrix
    derM = cat(3, derb0mu1M, derb1mu1M, derb0mu2M, derb1mu2M, derb0hM, derb1hM, derpM);

    %------Non zero Second Derivatives--------------

    %--entry (9,1)---
     %            UU         UM             MU           UH       HU       HM       MH       MM        HH         
    derp_b0mu1M  = [0         0              0           0        0        0        0        0         0;  %UU
                    0         0              0           0        0        0        0        0         0;  %UM
                    0         0              0           0        0        0        0        0         0;  %MU
                    0         0              0           1-mu2    0        0        mu2-1    0         0;  %UH
                    0         0              0           0        1-mu2    mu2-1    0        0         0;  %HU
                    0         0              0           0        0        0        0        0         0;  %HM
                    0         0              0           0        0        0        0        0         0;  %MH
                    0         0              0           0        0        0        0        0         0;  %MM
                    0         0              0           0        0        0        0        0         0;];%HH   

    %--entry (9,3)---            
    %             UU        UM       MU          UH       HU       HM       MH       MM        HH         
    derp_b0mu2M = [0         0        0           0        0        0        0        0         0;  %UU
                   0         0        0           0        0        0        0        0         0;  %UM
                   0         0        0           0        0        0        0        0         0;  %MU
                   0         0        0           -mu1     0        0        mu1      0         0;  %UH
                   0         0        0           0        -mu1     mu1      0        0         0;  %HU
                   0         0        0           0        0        0        0        0         0;  %HM
                   0         0        0           0        0        0        0        0         0;  %MH
                   0         0        0           0        0        0        0        0         0;  %MM
                   0         0        0           0        0        0        0        0         0;];%HH              

    %--entry (9,2)---            
    derp_b1mu1M = derp_b0mu1M * t;
    %--entry (9,4)---    
    derp_b1mu2M = derp_b0mu2M * t;


    %--entry (3,3)---   
    %                  UU        UM             MU          UH       HU       HM       MH       MM        HH         
    derb0mu2_b0mu2M = [2         -2             -2          0        0        0        0        2         0;  %UU
                       0         0              0           0        0        0        0        0         0;  %UM
                       0         0              0           0        0        0        0        0         0;  %MU
                       0         0              0           0        0        0        0        0         0;  %UH
                       0         0              0           0        0        0        0        0         0;  %HU
                       0         0              0           0        0        0        0        0         0;  %HM
                       0         0              0           0        0        0        0        0         0;  %MH
                       0         0              0           0        0        0        0        0         0;  %MM
                       0         0              0           0        0        0        0        0         0;];%HH   
    %--entry(3,4)---  
    derb0mu2_b1mu2M = derb0mu2_b0mu2M * t;
    %--entry(4,4)---
    derb1mu2_b1mu2M = derb0mu2_b1mu2M * t;


    %compress the second derivative matrices in one
    %considering the symmetries of the second derivatives
    secondDerM = zeros(9, 9, 7, 7);
    secondDerM(:,:,7,1) = derp_b0mu1M;
    secondDerM(:,:,1,7) = derp_b0mu1M;
    secondDerM(:,:,7,3) = derp_b0mu2M;
    secondDerM(:,:,3,7) = derp_b0mu2M;
    secondDerM(:,:,7,2) = derp_b1mu1M;
    secondDerM(:,:,2,7) = derp_b1mu1M;
    secondDerM(:,:,7,4) = derp_b1mu2M;
    secondDerM(:,:,4,7) = derp_b1mu2M;
    secondDerM(:,:,3,3) = derb0mu2_b0mu2M;
    secondDerM(:,:,4,4) = derb1mu2_b1mu2M;
    secondDerM(:,:,3,4) = derb0mu2_b1mu2M;
    secondDerM(:,:,4,3) = derb0mu2_b1mu2M;
end


end

