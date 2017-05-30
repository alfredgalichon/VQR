function [ Xv ] = vec( X)
%VEC Summary of this function goes here
%   Detailed explanation goes here

[r,q] = size(X);

Xv = reshape(X,r*q,1);
end

