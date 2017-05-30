function [pi,psi,b, val ] = MKQRTp( X,Y,U,mu,nu,solver )
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
c=-(vec(U*Y'))';
% c=-(kron(Y,U)*vec(eye(d)))';
A1=kron(sparse(1:n,1:n,1),sparse(ones(m,1)'));
A2=kron(X',sparse(1:m,1:m,1));
d1=vec(nu');
d2=vec(mu*xbar);
e=ones(m*n,1);
A=vertcat(A1,A2);
d=vertcat(d1,d2);
pi_init=vec(mu*nu');

% -------------------------------------------------------------------
%-------------------  solves the LP ---------------------------------
if solver=='gurobi'
    
    model.A=sparse(A);
    model.obj=c';
    model.modelsense='min';
    model.rhs=d;
    model.ub=e;
    model.sense='=';
    model.start=pi_init;
    result = gurobi(model);
    if strcmp(result.status, 'OPTIMAL')
        pivec=result.x;
        Lvec=result.pi';
    else
        fprintf('Optimization returned status: %s\n', result.status);
    end
else if solver=='knitro'
        [pivec, lambda, exitflag] = ktrlinklp (pi_init, c', [], [], sparse(A), d, zeros(m*n,1), e );
        Lvec=-lambda.eqlin';
        
    else
        
        warning('off','MATLAB:rankDeficientMatrix');
        opts1= optimset('display','off');
        %
        [pivec,Lvec,dum]=lp_fnm(A, c, d, e, pi_init);
        
        % Solve a linear program by the interior point method:
        % min(c * z), s.t. A * z = d and 0 < z < e
    end
end
% -------------------------------------------------------------------
% ------ reshape sols into matrix form and computes results  --------

pi=reshape(pivec,m,n);
L1vec=Lvec(1:n);
L2vec=Lvec((n+1):(n+m*r));
L1=reshape(L1vec,1,n);
L2=reshape(L2vec,m,r);

psi=-L1'; % beta is the regression coefficient
b=-L2;   % P is the positive part of Y-X.beta
val=trace(U'*pi*Y);           % value of the program
end
