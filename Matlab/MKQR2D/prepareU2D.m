function [ U_prov,m_prov,mu_prov ] = prepareU2D( T )
       [U1m,U2m]=meshgrid(T,T);
       U_prov=[vec(U1m) vec(U2m)];
       m_prov=size(U_prov,1);
       mu_prov=ones(m_prov,1)/m_prov;
end

