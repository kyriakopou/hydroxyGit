function [x0, finalIntLB, finalIntUB] = unifDependRnd(lbConst, ubConst, lastDay)

x0 = zeros(1, 7);
x0Const = unifrnd(lbConst, ubConst);
%put them in x0 array
x0(1) = x0Const(1);
x0(3) = x0Const(2);
x0(5) = x0Const(3);
x0(7) = x0Const(4);

lbGrad = -x0Const / lastDay;
ubGrad = (1 - x0Const) / lastDay;

finalIntLB = [lbConst(1) lbGrad(1) lbConst(2) lbGrad(2) lbConst(3) lbGrad(3) lbConst(4)];
finalIntUB = [ubConst(1) ubGrad(1) ubConst(2) ubGrad(2) ubConst(3) ubGrad(3) ubConst(4)];

%get the gradients of the efficiencies
x0Grad = unifrnd(lbGrad, ubGrad);
%x0 grad
x0(2) = x0Grad(1);
x0(4) = x0Grad(2);
x0(6) = x0Grad(3);





end