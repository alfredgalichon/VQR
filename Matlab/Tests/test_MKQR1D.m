[ret, compname] = system('hostname'); 
duma=double('AG-VAIO');
dumb=double(compname);
dumc=dumb(1:7);
homeComputer=max(duma==dumc);
if homeComputer
    solver='gurobi';
    path='C:\gurobi560\win64\matlab';
    step = 0.01;
else
    step = 0.001;
    solver='gurobi';
    path='C:\gurobi560\win32\matlab';

end
solver=initSolver(solver,path);
%step may go in practice from 0.2 to 0.001 (in d=1);

[X,Y,n,d,r]=LaunchData('Engel.xls');
if d~=1
        error('d has to be 1 or 2');
end

nu=ones(n,1)/n;
T=  (0:step:1)';
U=T;
m=size(U,1);
mu=ones(m,1)/m;
[pi,psi,b, val ] = MKQRTp( X,Y,U,mu,nu,solver );
beta=ComputeBeta1D( mu,b );

x_quantiles = quantile(X,0:0.25:1); % valeurs desirees
Nb_display= size(x_quantiles,1);
yhat_MKQR=x_quantiles*beta';

for i=1:Nb_display
    plot(T(2:m-1),yhat_MKQR(i,2:m-1),'Color',[0 0.2*i 0]);
    hold on;
end
hold off

