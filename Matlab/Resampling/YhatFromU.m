function [ys,j] = YhatFromU( X,beta,U, u0 )
[v,j0]=min((U(:,1)-u0(1)) .^ 2 + (U(:,2)-u0(2)).^ 2);
ytemp=yhat2D(X,beta);
ys=ytemp(:,j,:);
end

