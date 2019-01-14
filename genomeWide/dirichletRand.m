function r = dirichletRand(a, n)
%DIRICHLETRAND Summary of this function goes here
%  sample from dirichlet distribution n times
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

end

