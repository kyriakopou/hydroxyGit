function [A, derA, secondDerA] = demethylation(d, t, derCompFlag)
%This is active demethylation that models H->U 

%    UU      UM        MU         UH       HU       HM     MH     MM        HH  
A = [1       0         0          0        0        0      0      0         0            ;  %UU
     0       1         0          0        0        0      0      0         0            ;  %UM
     0       0         1          0        0        0      0      0         0            ;  %MU
     d       0         0          1-d      0        0      0      0         0            ;  %UH
     d       0         0          0        1-d      0      0      0         0            ;  %HU
     0       d         0          0        0        1-d    0      0         0            ;  %HM
     0       0         d          0        0        0      1-d    0         0            ;  %MH
     0       0         0          0        0        0      0      1         0            ;  %MM
     d^2     0         0          d*(1-d)  d*(1-d)  0      0      0         (1-d)^2      ;];%HH

derA = [];
secondDerA = [];
 
if derCompFlag
        %      UU      UM        MU         UH       HU       HM     MH     MM        HH  
    derb0dA = [0       0         0          0        0        0      0      0         0            ;  %UU
               0       0         0          0        0        0      0      0         0            ;  %UM
               0       0         0          0        0        0      0      0         0            ;  %MU
               1       0         0          -1       0        0      0      0         0            ;  %UH
               1       0         0          0        -1       0      0      0         0            ;  %HU
               0       1         0          0        0        -1     0      0         0            ;  %HM
               0       0         1          0        0        0      -1     0         0            ;  %MH
               0       0         0          0        0        0      0      0         0            ;  %MM
               2*d     0         0          1-2*d    1-2*d    0      0      0         -2+2*d       ;];%HH

    derb1dA = derb0dA * t;   

    derb0mu1A = zeros(9,9);
    derb1mu1A = zeros(9,9);
    derb0mu2A = zeros(9,9);
    derb1mu2A = zeros(9,9);
    derb0hA = zeros(9,9); 
    derb1hA = zeros(9,9);      
    derpA = zeros(9,9); 

    %concatenate  all derivative matrices into a 3 dimensional matrix
    derA = cat(3, derb0mu1A, derb1mu1A, derb0mu2A, derb1mu2A, derb0hA, derb1hA, derb0dA, derb1dA, derpA);


    %------Non zero Second Derivatives--------------
    %entry (7,7)
            %      UU      UM        MU       UH       HU       HM     MH     MM        HH  
    derb0d_b0dA = [0       0         0        0        0        0      0      0         0            ;  %UU
                   0       0         0        0        0        0      0      0         0            ;  %UM
                   0       0         0        0        0        0      0      0         0            ;  %MU
                   0       0         0        0        0        0      0      0         0            ;  %UH
                   0       0         0        0        0        0      0      0         0            ;  %HU
                   0       0         0        0        0        0      0      0         0            ;  %HM
                   0       0         0        0        0        0      0      0         0            ;  %MH
                   0       0         0        0        0        0      0      0         0            ;  %MM
                   2       0         0        -2       -2       0      0      0         2            ;];%HH

    %entry(7,8)
    derb0d_b1dA = derb0d_b0dA * t;
    %entry(8,8)
    derb1d_b1dA = derb0d_b1dA * t;

    %compress the second derivative matrices in one
    %considering the symmetries of the second derivatives
    secondDerA = zeros(9, 9, 9, 9);
    secondDerA(:,:,7,7) = derb0d_b0dA;
    secondDerA(:,:,7,8) = derb0d_b1dA;
    secondDerA(:,:,8,7) = derb0d_b1dA;
    secondDerA(:,:,8,8) = derb1d_b1dA;
end

end

