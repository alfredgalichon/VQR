function [ beta,U,m,mu,pi,b ] = ComputeBetaEtAl2D( b_prov,T,U_prov,pi_prov,step )
U1_prov=U_prov(:,1);
U2_prov=U_prov(:,2);
[m_prov,r]=size(b_prov);
[m_prov,n]=size(pi_prov);
nonzind=find(U1_prov~=0 & U2_prov~=0);
U1=U1_prov(nonzind);
U2=U2_prov(nonzind);
U=[U1 U2];
m=size(U,1);
mu=ones(m,1)/m;
beta=zeros(m,r,2);
for k=1:r
    [D1bk_prov,D2bk_prov]=grad2D(b_prov(:,k),T,U1_prov, U2_prov,step);
    beta(:,k,1)=D1bk_prov(nonzind);
    beta(:,k,2)=D2bk_prov(nonzind);
    
end
pi=zeros(m,n);
for i=1:n
    pi(:,i)=pi_prov(nonzind,i);
end
b=b_prov(nonzind);
end

