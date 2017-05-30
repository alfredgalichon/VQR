save=true;
[ret, compname] = system('hostname'); 
duma=double('AG-VAIO');
dumb=double('PST5370');
dum=double(compname);
dumc=dum(1:7);
homeComputer=max(duma==dumc);
officeComputer=max(dumb==dumc);
if homeComputer
    solver='gurobi';
 	path='C:\gurobi560\win64\matlab';
else if officeComputer
    solver='gurobi';
    path='C:\gurobi560\win32\matlab';
    end
end

solver=initSolver(solver,path);
%step may go in practice from 0.2 to 0.001 (in d=1);

[X,Yprov,n,d,r]=LaunchData('MVEngel.xls');

% -----------------------------
% Phase 1: two univariate mu-QR
% -----------------------------
Y=[Yprov(:,1),(Yprov(:,3)+Yprov(:,4))];
step = 0.01;
nu=ones(n,1)/n;
T=  (0:step:1)';
lT=size(T,1);
U=T;
m=size(U,1);
mu=ones(m,1)/m;


% -----------------------------
% Phase 2: one bivariate mu-QR
% -----------------------------

step = 0.1; 
nu=ones(n,1)/n;
T=  (0:step:1)';
lT=size(T,1);
U=T;
m=size(U,1);
mu=ones(m,1)/m;

[ U_prov,m_prov,mu_prov ] = prepareU2D( T );
[pi_prov,psi_prov,b_prov, val ] = MKQRTp( X,Y,U_prov,mu_prov,nu,solver );

[pi_prov_test,psi_prov_test,b_prov_test, val_test ] = MKQRTd( X,Y,U_prov,mu_prov,nu,solver );
[ beta,U,m,mu,pi,b ] = ComputeBetaEtAl2D( b_prov,T,U_prov,pi_prov,step );

