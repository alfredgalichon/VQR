function [ys] = yhat2D( xs, beta )
[m,r,d]=size(beta);
[p,rbis]=size(xs);
if or(r~=rbis,d~=2)
    error('Error in yhat2D');
else
    ys=zeros(p,m,d);
    for i=1:p
        for j=1:m
        for l=1:d
            ys(i,j,l)=beta(j,:,l)*xs(i,:)';
            
        end
        end
    end
    
end

end

