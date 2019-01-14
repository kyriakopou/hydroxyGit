function [F,dF]=dexpme(M,dM)
% Eigenvalues/eigenvectors decomposition - usage of Matlab eig
[U,Q]=eig(M);
eQr=exp(diag(Q));
H=U\dM*U;
P=diag(diag(H).*eQr);
for i=1:length(H)-1
    for j=i+1:length(H)
        RQ=(eQr(i)-eQr(j))/(Q(i,i)-Q(j,j));
        P(i,j)=H(i,j)*RQ; P(j,i)=H(j,i)*RQ;
    end
end
% by Lubomír Bran?ík, 2008
F=U*diag(eQr)/U;
dF=U*P/U;