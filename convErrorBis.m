function Cbis = convErrorBis(errorsBis)

%assign the error probabilities
c = 1-errorsBis(1);
d = 1-errorsBis(2);
e = 1-errorsBis(3);

%       TT           TC           CT           CC     
Cbis = [c^2          c*(1-c)      c*(1-c)      (1-c)^2        ;  %UU
        c*(1-d)      c*d          (1-c)*(1-d)  d*(1-c)        ;  %UM
        (1-d)*c      (1-c)*(1-d)  d*c          d*(1-c)        ;  %MU
        c*(1-e)      c*e          (1-c)*(1-e)  (1-c)*e        ;  %UH
        (1-e)*c      (1-e)*(1-c)  c*e          e*(1-c)        ;  %HU
        (1-e)*(1-d)  (1-e)*d      e*(1-d)      e*d            ;  %HM
        (1-e)*(1-d)  (1-d)*e      (1-e)*d      e*d            ;  %MH
        (1-d)^2      d*(1-d)      d*(1-d)      d^2            ;  %MM
        (1-e)^2      (1-e)*e      e*(1-e)      e^2           ];  %HH
   


end