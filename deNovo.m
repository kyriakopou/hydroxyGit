function [Md, derMd, secondDerMd] = deNovo(mu2, t, derCompFlag)
%de-Novo Methylation Process Matrix
%mu2 is de-Novo methylation probability


tot = 2*mu2-mu2^2;   

%     UU        UM             MU          UH       HU       HM       MH        MM        HH                   
Md = [1-tot     mu2-mu2^2      mu2-mu2^2   0        0        0        0         mu2^2     0;  %UU
      0         1-mu2          0           0        0        0        0         mu2       0;  %UM
      0         0              1-mu2       0        0        0        0         mu2       0;  %MU
      0         0              0           1-mu2    0        0        mu2       0         0;  %UH
      0         0              0           0        1-mu2    mu2      0         0         0;  %HU
      0         0              0           0        0        1        0         0         0;  %HM
      0         0              0           0        0        0        1         0         0;  %MH
      0         0              0           0        0        0        0         1         0;  %MM
      0         0              0           0        0        0        0         0         1]; %HH 
derMd = [];
secondDerMd = [];
  
if derCompFlag               
    %             UU        UM             MU          UH       HU       HM       MH       MM        HH         
    derb0mu2Md = [-2+2*mu2  1-mu2*2        1-mu2*2     0        0        0        0        mu2*2     0;  %UU
                  0         -1             0           0        0        0        0        1         0;  %UM
                  0         0              -1          0        0        0        0        1         0;  %MU
                  0         0              0           -1       0        0        1        0         0;  %UH
                  0         0              0           0        -1       1        0        0         0;  %HU
                  0         0              0           0        0        0        0        0         0;  %HM
                  0         0              0           0        0        0        0        0         0;  %MH
                  0         0              0           0        0        0        0        0         0;  %MM
                  0         0              0           0        0        0        0        0         0]; %HH   


    derb1mu2Md = derb0mu2Md * t;
    derb0mu1Md = zeros(9, 9); 
    derb1mu1Md = derb0mu1Md * t;
    derb0hMd = zeros(9, 9);
    derb1hMd = derb0hMd * t;
    derb0dMd = zeros(9, 9);
    derb1dMd = derb0dMd * t;     
    derpMd = zeros(9, 9); 



    %concatenate  all derivative matrices into a 3 dimensional matrix
    derMd = cat(3, derb0mu1Md, derb1mu1Md, derb0mu2Md, derb1mu2Md, derb0hMd, derb1hMd, derb0dMd, derb1dMd, derpMd);   


    %------Non zero Second Derivatives--------------
    %--entry (3,3)---
    %                   UU        UM             MU          UH       HU       HM       MH       MM        HH         
    derb0mu2_b0mu2Md = [2         -2             -2          0        0        0        0        2         0;  %UU
                        0         0              0           0        0        0        0        0         0;  %UM
                        0         0              0           0        0        0        0        0         0;  %MU
                        0         0              0           0        0        0        0        0         0;  %UH
                        0         0              0           0        0        0        0        0         0;  %HU
                        0         0              0           0        0        0        0        0         0;  %HM
                        0         0              0           0        0        0        0        0         0;  %MH
                        0         0              0           0        0        0        0        0         0;  %MM
                        0         0              0           0        0        0        0        0         0]; %HH  

    %--entry (4,3)---
    derb1mu2_b0mu2Md = derb0mu2_b0mu2Md * t;
    %--entry (4,4)---
    derb1mu2_b1mu2Md = derb1mu2_b0mu2Md * t;

    %compress the second derivative matrices in one
    %considering the symmetries of the second derivatives
    secondDerMd = zeros(9, 9, 9, 9);
    secondDerMd(:,:,3,3) = derb0mu2_b0mu2Md;
    secondDerMd(:,:,4,3) = derb1mu2_b0mu2Md;
    secondDerMd(:,:,3,4) = derb1mu2_b0mu2Md;
    secondDerMd(:,:,4,4) = derb1mu2_b1mu2Md;
end

end