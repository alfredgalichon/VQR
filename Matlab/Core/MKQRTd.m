function [pi,psi,b, val ] = MKQRTd( X,Y,U,mu,nu,solver )
% Monge-Kantorovich (transportation version)
% computes standard classical regression
% programmed by Alfred Galichon, with research assistance from Jinxin He
% reference: Roger Koenker; Gilbert Bassett, Jr. (1978). "Regression
% Quantiles". Econometrica, Vol. 46, No. 1. (Jan., 1978), pp. 33-50.

% --------------- construction of the constant matrices --------------
LARGE=1000000;
[n,d]=size(Y);
r=size(X,2);
[m, dprime]=size(U);
if d ~= dprime
    error('Error: dimensions of Y and U must agree');
    return
end

xbar=nu'*X;
% -------------------------------------------------------------------
% --------------- vectorization of the matrices for LP --------------
c1=kron(1,nu');
c2=kron(xbar,mu');
c=horzcat(c1,c2);
A1=kron(sparse(ones(m,1)),sparse(1:n,1:n,1));
A2=kron(sparse(1:m,1:m,1),X)*transv(m,r);
d=vec(Y*U');
A=horzcat(A1,A2);

psi_init=max(U*Y')';
bv_init=zeros(m*r,1);
z_init=vertcat(psi_init,bv_init);
% -------------------------------------------------------------------
%-------------------  solves the LP ---------------------------------
if solver=='gurobi'
    
    model.A=sparse(A);
    model.obj=c';
    model.modelsense='min';
    model.rhs=d;
    model.sense='>';
    model.start=z_init;
    result = gurobi(model);
    if strcmp(result.status, 'OPTIMAL')
        L1=result.x;
        pivec=result.pi;
    else
        fprintf('Optimization returned status: %s\n', result.status);
    end
else if solver=='knitro'
        ktrlinklp (pi_init, c', [], [], sparse(A), d, zeros(m*n,1), e );
        0/0;
        
    else
        
        warning('off','MATLAB:rankDeficientMatrix');
        opts1= optimset('display','off');
        %
        0/0;
        [pivec,Lvec,dum]=lp_fnm(A, c, d, e, pi_init);
        
        [pivec_bis, lambda, exitflag] = ktrlinklp (pi_init, c', [], [], sparse(A), d, zeros(m*n,1), e );
        Lvec_bis=-lambda.eqlin';
        probleme ici avec Lvec_bus neq Lvec
        % Solve a linear program by the interior point method:
        % min(c * z), s.t. A * z = d and 0 < z < e
    end
end
% -------------------------------------------------------------------
% ------ reshape sols into matrix form and computes results  --------

pi=reshape(pivec,n,m)';
L1vec=L1(1:n);
L2vec=L1((n+1):(n+m*r));
L1=reshape(L1vec,1,n);
L2=reshape(L2vec,m,r);

psi=L1'; % beta is the regression coefficient
b=L2;   % P is the positive part of Y-X.beta
val=trace(U'*pi*Y)           % value of the program

% -------------------------------------------------------------------
% ------ second round to guarantee uniqueness  ----------------------

cbis=horzcat(c1,-c2);


A1=kron(sparse(ones(m,1)),sparse(1:n,1:n,1));
A2=kron(sparse(1:m,1:m,1),X)*transv(m,r);
Apos= eye(n,(n+m*r));
dpos=zeros(n,1);
Abis=vertcat(A,Apos,-c);
dbis=vertcat(d,dpos,-val);

psi_init=L1vec;
bv_init=L2vec;
z_init=vertcat(psi_init,bv_init);


% -------------------------------------------------------------------
%-------------------  solves the LP again ---------------------------
if solver=='gurobi'
    
    model.A=sparse(Abis);
    model.obj=cbis';
    model.modelsense='min';
    model.rhs=dbis;
    model.sense='>';
    model.start=z_init;
    result = gurobi(model);
    if strcmp(result.status, 'OPTIMAL')
        L1=result.x;
        pivec=result.pi;
    else
        fprintf('Optimization returned status: %s\n', result.status);
    end
else probleme ici;
end
% -------------------------------------------------------------------
% ------ reshape sols into matrix form and computes results  --------

pi=reshape(pivec(1:n*m),n,m)';
L1vec=L1(1:n);
L2vec=L1((n+1):(n+m*r));
L1=reshape(L1vec,1,n);
L2=reshape(L2vec,m,r);

psi=L1'; % beta is the regression coefficient
b=L2;   % P is the positive part of Y-X.beta
val=trace(U'*pi*Y)           % value of the program


end
