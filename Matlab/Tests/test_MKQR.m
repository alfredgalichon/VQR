solver='gurobi';
solver=initSolver(solver);
step = 0.01;
%step may go in practice from 0.2 to 0.001 (in d=1);

[X,Y,n,d,r]=LaunchData('GDPDebtDeficit.xls');
nu=ones(n,1)/n;
T=  (0:step:1)';
lT=size(T,1);

if d==1
    U=T;
    m=size(U,1);
    mu=ones(m,1)/m;
    [pi,psi,b, val ] = MKQRTp( X,Y,U,mu,nu,solver );
    beta=ComputeBeta1D( mu,b );
    
    x_quantiles = quantile(X,0:0.25:1); % valeurs desirees
    Nb_display= size(x_quantiles,1);
    yhat_MKQR=x_quantiles*beta';
else
    if d==2
        [ U_prov,m_prov,mu_prov ] = prepareU2D( T );
        [pi_prov,psi_prov,b_prov, val ] = MKQRTp( X,Y,U_prov,mu_prov,nu,solver );
        [ beta,U,m,mu,pi,b ] = ComputeBetaEtAl2D( b_prov,T,U_prov,pi_prov,step );

        xeval=[1.0000 2.33];
        PlotFixedx2D( xeval,beta,lT );
        [ Yhat,pihat ] = bootstrap( X,beta,pi );
        
    else
        error('d has to be 1 or 2');
    end
end

