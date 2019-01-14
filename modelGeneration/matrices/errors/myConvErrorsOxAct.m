function E = myConvErrorsOxAct(convProbs)
%MYCONVERRORSOXACT
%    E = MYCONVERRORSOXACT(C,D,F,G)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    16-Mar-2017 16:15:16

c = convProbs(1);
d = convProbs(2);
f = convProbs(3);
g = convProbs(4);

t2 = c-1.0;
t3 = d-1.0;
t4 = f-1.0;
t5 = g-1.0;
t6 = t2.*t3;
t7 = c.*d;
t8 = c.*f;
t9 = t2.*t4;
t10 = d.*f;
t11 = t3.*t4;
t12 = c.*g;
t13 = t2.*t5;
t14 = d.*g;
t15 = t3.*t5;
t16 = f.*g;
t17 = t4.*t5;
E = reshape([c.^2,-c.*t3,t8,t12,-c.*t3,t3.^2,-f.*t3,-g.*t3,t8,-f.*t3,f.^2,t16,t12,-g.*t3,t16,g.^2,-c.*t2,t7,-c.*t4,-c.*t5,t6,-d.*t3,t11,t15,-f.*t2,t10,-f.*t4,-f.*t5,-g.*t2,t14,-g.*t4,-g.*t5,-c.*t2,t6,-f.*t2,-g.*t2,t7,-d.*t3,t10,t14,-c.*t4,t11,-f.*t4,-g.*t4,-c.*t5,t15,-f.*t5,-g.*t5,t2.^2,-d.*t2,t9,t13,-d.*t2,d.^2,-d.*t4,-d.*t5,t9,-d.*t4,t4.^2,t17,t13,-d.*t5,t17,t5.^2],[16,4]);