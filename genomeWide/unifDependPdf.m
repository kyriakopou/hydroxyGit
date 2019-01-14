function pdf = unifDependPdf(x, lastDay)

xConst = [x(1) x(3) x(5)];
xGrad = [x(2) x(4) x(6)];
lbGrad = -xConst / lastDay;
ubGrad = (1 - xConst) / lastDay;

pdf = prod(unifpdf(xGrad, lbGrad, ubGrad));

end