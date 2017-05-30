save=true;
solver='gurobi';
path='C:\gurobi604\win64\matlab';

solver=initSolver(solver,path);
%step may go in practice from 0.2 to 0.001 (in d=1);

[X,Y,n,d,r]=LaunchData('MVEngel');

% -----------------------------
% Phase 1: two univariate mu-QR
% -----------------------------
step = 0.01;
nu=ones(n,1)/n;
T=  (0:step:1)';
lT=size(T,1);
U=T;
m=size(U,1);
mu=ones(m,1)/m;
for theDim=1:2
    [pi,psi,b, val ] = MKQRTd( X,Y(:,theDim),U,mu,nu,solver );
    beta=ComputeBeta1D( mu,b );
    
    x_quantiles = quantile(X,0.5); 
    %x_quantiles = quantile(X,0:0.25:1); % valeurs desirees
    Nb_display= size(x_quantiles,1);
    yhat_MKQR=x_quantiles*beta';
    
    for i=1:Nb_display
        plot(T(2:m-1),yhat_MKQR(i,2:m-1),'Color',[0 0.2*i 0]);
        hold on;
    end
    xlabel('u \in [0,1]')
    ylabel(strcat('y_',num2str(theDim)))
    title(strcat('$u \to x \cdot \beta^{QR} (u)$ for $Y_',num2str(theDim),'$'), 'Interpreter','latex')
    if save
        h=gcf;
        print(h,'-dpdf',strcat('Y',num2str(theDim),'.pdf')) 
    end
    hold off
    if theDim<2 
        figure
    end
end


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

lambda = 1;

[ U_prov,m_prov,mu_prov ] = prepareU2D( T );
[pi_prov,psi_prov,b_prov, val ] = MKQRTd( X,Y .* (ones(n,1) * [1, lambda]),U_prov,mu_prov,nu,solver );

[ beta,U,m,mu,pi,b ] = ComputeBetaEtAl2D( b_prov,T,U_prov,pi_prov,step );
beta(:,:,2)=beta(:,:,2)*lambda;

xmedian=quantile(X,0.5);
PlotFixedx2D( xmedian,beta,lT,save,'xbeta',' \beta_1(u_1, u_2)^\top x','', '}$, for x=(1, 883.99)' );

%xeval=[1 0]; % not really needed as takes the constant part of the
%regressor
%PlotFixedx2D( xeval,beta,lT,save,'beta',' \beta','1','}(u_1,u_2)$' );
%xeval=[0 1];
%PlotFixedx2D( xeval,beta,lT, save,'beta',' \beta','2','}(u_1,u_2)$');

%xmean=mean(X);
%PlotFixedx2D( xmean,beta,lT,save,'xbarbeta',' E[Y', '','}|U=u]$' );

%[ Yhat,pihat ] = bootstrap( X,beta,pi );
        
   
