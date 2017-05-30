function [ Yhat,pihat ] = bootstrap( X,beta,pif )
Yhatf=yhat2D( X, beta );
[n,m,d]=size(Yhatf);
[mbis,nbis]=size(pif);
if or(nbis~=n,mbis~=m)
    error('Error in bootstrap');
else
    Yhat=zeros(n*m,d);
    for l=1:d
        Yhat(:,l)=reshape(Yhatf(:,:,l),m*n,1);
    end
end
pihat=reshape(pif,m*n,1);

end

