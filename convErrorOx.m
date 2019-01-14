function COxbis = convErrorOx(errorsOx)

%assign the error probabilities
c = 1-errorsOx(1);
d = 1-errorsOx(2);
f = 1-errorsOx(3);

%         TT       TC           CT             CC     
COxbis = [c^2      c*(1-c)      c*(1-c)        (1-c)^2        ;  %UU
          c*(1-d)  c*d          (1-c)*(1-d)    d*(1-c)        ;  %UM
          (1-d)*c  (1-c)*(1-d)  d*c            d*(1-c)        ;  %MU
          c*f      c*(1-f)      (1-c)*f        (1-c)*(1-f)    ;  %UH
          c*f      (1-c)*f      c*(1-f)        (1-c)*(1-f)    ;  %HU
          f*(1-d)  f*d          (1-f)*(1-d)    (1-f)*d        ;  %HM
          f*(1-d)  (1-f)*(1-d)  f*d            (1-f)*d        ;  %MH
          (1-d)^2  (1-d)*d      d*(1-d)        d^2            ;  %MM
          f^2      f*(1-f)      f*(1-f)        (1-f)^2       ];  %HH
   


end