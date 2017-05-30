function [X,Y,n,d,r]=LaunchData(datasetName)
Y = csvread(strcat(datasetName,'Y.csv'));
[n,d]=size(Y);
X = [ones(n,1), csvread(strcat(datasetName,'X.csv'))];
[ndum,r]=(size(X));
end