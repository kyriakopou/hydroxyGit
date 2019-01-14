function [H, derH, secondDerH] = hydroxylation(h, t, derCompFlag)
%HydroxYlation Process Matrix

%     UU        UM        MU        UH        HU      HM           MH           MM           HH
H =  [1         0         0         0         0       0            0            0            0     ;     %UU
      0         1-h       0         h         0       0            0            0            0     ;     %UM
      0         0         1-h       0         h       0            0            0            0     ;     %MU
      0         0         0         1         0       0            0            0            0     ;     %UH
      0         0         0         0         1       0            0            0            0     ;     %HU
      0         0         0         0         0       1-h          0            0            h     ;     %HM
      0         0         0         0         0       0            (1-h)        0            h     ;     %MH
      0         0         0         0         0       h-h^2        h-h^2        (1-h)*(1-h)  h^2   ;     %MM
      0         0         0         0         0       0            0            0            1     ;];   %HH  
derH = []; 
secondDerH = [];

if derCompFlag
        %     UU       UM        MU        UH        HU      HM           MH            MM        HH
    derb0hH = [0         0         0         0         0       0            0             0         0     ;     %UU
               0         -1        0         1         0       0            0             0         0     ;     %UM
               0         0         -1        0         1       0            0             0         0     ;     %MU
               0         0         0         0         0       0            0             0         0     ;     %UH
               0         0         0         0         0       0            0             0         0     ;     %HU
               0         0         0         0         0       -1           0             0         1     ;     %HM
               0         0         0         0         0       0            -1            0         1     ;     %MH
               0         0         0         0         0       1-2*h        1-2*h         -2+2*h    2*h   ;     %MM
               0         0         0         0         0       0            0             0         0     ;];   %HH


    derb1hH = derb0hH * t;

    derb0mu1H = zeros(9,9);
    derb1mu1H = zeros(9,9);
    derb0mu2H = zeros(9,9);
    derb1mu2H = zeros(9,9);
    derpH = zeros(9,9);    
    %concatenate  all derivative matrices into a 3 dimensional matrix
    derH = cat(3, derb0mu1H, derb1mu1H, derb0mu2H, derb1mu2H, derb0hH, derb1hH, derpH);


    %------Non zero Second Derivatives--------------

    %entry (5,5) 
       %     UU       UM        MU        UH        HU      HM           MH            MM        HH
    derb0h_b0hH = [0         0         0         0         0       0            0             0         0     ;     %UU
                   0         0         0         0         0       0            0             0         0     ;     %UM
                   0         0         0         0         0       0            0             0         0     ;     %MU
                   0         0         0         0         0       0            0             0         0     ;     %UH
                   0         0         0         0         0       0            0             0         0     ;     %HU
                   0         0         0         0         0       0            0             0         0     ;     %HM
                   0         0         0         0         0       0            0             0         0     ;     %MH
                   0         0         0         0         0       -2           -2            2         2     ;     %MM
                   0         0         0         0         0       0            0             0         0     ;];   %HH

    %entry (5,6)
    derb0h_b1hH = derb0h_b0hH * t;
    derb1h_b1hH = derb0h_b1hH * t;

    %compress the second derivative matrices in one 4 dimensional matrix
    %considering the symmetries of the second derivatives
    secondDerH = zeros(9, 9, 7, 7);
    secondDerH(:,:,5,5) = derb0h_b0hH;
    secondDerH(:,:,5,6) = derb0h_b1hH;
    secondDerH(:,:,6,5) = derb0h_b1hH;
    secondDerH(:,:,6,6) = derb1h_b1hH;
end

end

