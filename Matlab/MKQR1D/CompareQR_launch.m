function []=CompareQR_launch(file);
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
EPS=0.0001;
%step = 0.2;
%step = 0.1;
%step = 0.05;
%  step = 0.01;
% step = 0.005;
step = 0.01;
%step = 0.0005;
T=  EPS:step:(1+EPS-step);
m=length(T);
T=T';


D = eye(m);
for i=1:m-1 D(i+1,i)=-1; end
U=vertcat(T,1);
U(1,:)=[];
mu=D*U;

S=ones(m);
for i=1:m
    for j=(i+1):m
        S(i,j)=0;
    end
end

[X,Y,n,d,r]=LaunchData(file);
nu=ones(n,1)/n;

x_quantiles = quantile(X,0:0.25:1); % valeurs desirees
Nb_display= size(x_quantiles,1);
[ Z_QR,P_QR,beta_QR, val_QR ]= QRp( X,Y,T,nu, solver );
yhat_QR=x_quantiles*beta_QR';
[pi_MKQR,psi_MKQR,b_MKQR, val_MQKR ] = MKQRTp( X,Y,U,mu,nu , solver);
beta_MKQR=inv(diag(mu))*D*b_MKQR;
yhat_MKQR=x_quantiles*beta_MKQR';

%[ Z_KNQR,Gamma_KNQR, P_KNQR,beta_KNQR, val_KNQR ]= KNQRp( X,Y,T,nu,solver );
%yhat_KNQR=x_quantiles*beta_KNQR';

for i=1:Nb_display
    plot(T(2:m-1),yhat_QR(i,2:m-1),'Color',[0.2*i 0 0]);
    hold on;
    plot(T(2:m-1),yhat_MKQR(i,2:m-1),'Color',[0 0.2*i 0]);
    hold on;
    %plot(T(2:m-1),yhat_KNQR(i,2:m-1),'Color',[0.2*i 0 0]);
    %hold on;
end
hold off

end