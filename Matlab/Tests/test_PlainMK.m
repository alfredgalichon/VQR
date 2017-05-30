solver='gurobi'
solver=initSolver(solver);
step = 0.01;

[X,Y,n,d,r]=LaunchData('GDPDebtDeficitNoReg.xls');
nu=ones(n,1)/n;
T=  (0:step:1)';
lT=size(T,1);

if d~=2
    error('Not the right dimension of the data')
else
        [ U_prov,m_prov,mu_prov ] = prepareU2D( T );
        [pi_prov,psi_prov,b_prov, val ] = MKQRTp( X,Y,U_prov,mu_prov,nu,solver );
        [ beta,U,m,mu,pi,b ] = ComputeBetaEtAl2D( b_prov,T,U_prov,pi_prov,step );
        surf(reshape(b,lT-1,lT-1));
        xeval=[1.0000 ];
        PlotFixedx2D( xeval,beta,lT );
        [ Yhat,pihat ] = bootstrap( X,beta,pi );
end

