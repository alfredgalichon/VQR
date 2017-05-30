function [ Z,P,beta, val ] = QRp( X,Y,T,nu, solver )
% Classical Quantile Regression
% computes standard classical regression 
% programmed by Alfred Galichon, with research assistance from Jinxin He
% reference: Roger Koenker; Gilbert Bassett, Jr. (1978). "Regression
% Quantiles". Econometrica, Vol. 46, No. 1. (Jan., 1978), pp. 33-50.

% --------------- construction of the constant matrices --------------
[n,d]=size(Y);
r=size(X,2);
m=length(T);
U=vertcat(T,1);
U(1,:)=[];
D = eye(m);
for i=1:m-1 D(i+1,i)=-1; end
mu=D*U;
xbar=nu'*X;
% -------------------------------------------------------------------
% --------------- vectorization of the matrices for LP --------------
c=-kron(Y,mu);
A=kron(X',sparse(1:m,1:m,1));
d=vec((ones(m,1)-T)*xbar);
e=vec(ones(m,1)*nu');
z_init=vec((ones(m,1)-T)*nu');
% -------------------------------------------------------------------
%-------------------  solves the LP ---------------------------------
if solver=='gurobi'

    model.A=sparse(A);
    model.obj=c;
    model.modelsense='min';
    model.rhs=d;
    model.ub=e;
    model.sense='=';
    model.start=z_init;
    result = gurobi(model);
    if strcmp(result.status, 'OPTIMAL')
        Zvec=result.x;
        L1vec=result.pi';
        L2vec=-(Zvec>0)' .* result.rc'; 
    else
        fprintf('Optimization returned status: %s\n', result.status);
    end
else if solver=='knitro'
        [Zvec, lambda, exitflag] = ktrlinklp (z_init, c, [], [], sparse(A), d, zeros(m,n), e );
        L1vec=-lambda.eqlin';
        L2vec=lambda.upper';
        
    else
        warning('off','MATLAB:rankDeficientMatrix');
        opts1= optimset('display','off');
        [Zvec,L1vec,L2vec]=lp_fnm(A, c', d, e, z_init);        
        % Solves a linear program by the interior point method:
        % min(c * z), s.t. A * z = d and 0 < z < e
    end
end
% -------------------------------------------------------------------
% ------ reshape sols into matrix form and computes results  --------

Z=reshape(Zvec,m,n);
L1=reshape(L1vec,m,r); % Lagrange mult associated to equality constraint
L2=reshape(L2vec,m,n); % Lagrange mult associated to inequality constraint

beta=-inv(diag(mu))*L1; % beta is the regression coefficient
P=-L2'*inv(diag(mu));   % P is the positive part of Y-X.beta
val=mu'*Z*Y;           % value of the program
end
