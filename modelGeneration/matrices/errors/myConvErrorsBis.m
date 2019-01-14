function E = myConvErrorsBis(convProbs)
%MYCONVERRORSBIS
%    E = MYCONVERRORSBIS(C,D,E)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    13-Mar-2017 17:16:08

c = convProbs(1);
d = convProbs(2);
e = convProbs(3);

t2 = c-1.0;
t3 = d-1.0;
t4 = e-1.0;
t5 = t2.*t3;
t6 = c.*d;
t7 = t2.*t4;
t8 = c.*e;
t9 = t3.*t4;
t10 = d.*e;
E = reshape([c.^2,-c.*t3,-c.*t4,-c.*t3,t3.^2,t9,-c.*t4,t9,t4.^2,-c.*t2,t6,t8,t5,-d.*t3,-e.*t3,t7,-d.*t4,-e.*t4,-c.*t2,t5,t7,t6,-d.*t3,-d.*t4,t8,-e.*t3,-e.*t4,t2.^2,-d.*t2,-e.*t2,-d.*t2,d.^2,t10,-e.*t2,t10,e.^2],[9,4]);
