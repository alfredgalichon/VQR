function [ T ] = transv( m,n)
%VEC Summary of this function goes here
%   Detailed explanation goes here

k=1:(m*n);
l=arrayfun(@transf,k);
s=ones(1,m*n);
T=sparse(k,l,s);
    


function [q]=transf(p)
  i=mod(p-1,n)+1;
  j=floor((p-1)/n)+1;
  q=m*(i-1)+j;
end


end


