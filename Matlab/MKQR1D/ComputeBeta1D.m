function [ beta ] = ComputeBeta1D( mu,b )
[m,dum]=size(mu);
D = eye(m);
for i=1:m-1 D(i+1,i)=-1; end
 beta=inv(diag(mu))*D*b;
end

