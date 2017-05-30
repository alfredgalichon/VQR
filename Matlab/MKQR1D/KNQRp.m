function [ Z,Gamma, P,beta, val ] = KNQRp( X,Y,T,nu ,solver)
% Koenker and Ng's constrained Quantile Regression
% computes standard classical regression 
% programmed by Alfred Galichon, with research assistance from Jinxin He
% reference: Koenker, R., P. Ng (2005), "Inequality constrained quantile 
% regression." Sankhy? 67, no. 2, 418--440.

% --------------- construction of the constant matrices --------------
[n,d]=size(Y);
r=size(X,2);
m=length(T);
U=vertcat(T,1);
U(1,:)=[];
D = eye(m);
for i=1:m-1 D(i+1,i)=-1; end
D0=D;
D0(:,m)=[];
mu=D*U;
xbar=nu'*X;
LARGE=1000000000;
% --------------- constants specific to KNQR ------------------------
epsilon=0.5*min(mu(1)*T(1),mu(m)*T(m));
delta=zeros(m,1);
delta(1)=epsilon/mu(1);
delta(m)=epsilon/mu(m);

% -------------------------------------------------------------------
% --------------- vectorization of the matrices for LP --------------
c1=-kron(Y',mu');
c2=zeros(1,(m-1)*n);
c=horzcat(c1,c2);
A1=kron(X',eye(m));
A2=-kron(X',inv(diag(mu))*D0);
A=horzcat(A1,A2);
d=vec((ones(m,1)-T)*xbar);
e1=vec(ones(m,1)*nu');
e2=vec(LARGE*ones(m-1,n));
e=vertcat(e1,e2);
z1_init=vec((ones(m,1)+delta-T)*nu');
z2_init=vec(epsilon*ones(m-1,1)*nu');
z_init=vertcat(z1_init,z2_init);
% -------------------------------------------------------------------
%-------------------  solves the LP ---------------------------------
if solver=='gurobi'

    model.A=sparse(A);
    model.obj=c';
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
else
warning('off','MATLAB:rankDeficientMatrix');
opts1= optimset('display','off');
% 
[Zvec,L1vec,L2vec]=lp_fnm(A, c, d, e, z_init);
% Solve a linear program by the interior point method:
% min(c * z), s.t. A * z = d and 0 < z < e 
end
% -------------------------------------------------------------------
% ------ reshape sols into matrix form and computes results  --------

Z1vec=Zvec(1:m*n);
Z2vec=Zvec(m*n+1:(m*n+(m-1)*n));
Z=reshape(Z1vec,m,n);
Gamma=reshape(Z2vec,m-1,n)';
L2avec=L2vec(1:m*n);
L1=reshape(L1vec,m,r); % Lagrange mult associated to equality constraint
L2=reshape(L2avec,m,n); % Lagrange mult associated to inequality constraint

beta=-inv(diag(mu))*L1; % beta is the regression coefficient
P=-L2'*inv(diag(mu));   % P is the positive part of Y-X.beta
val=mu'*Z*Y;           % value of the program
end
