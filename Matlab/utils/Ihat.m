function [valf]=Ihat(u0,vect,U,sigma)
% kernel smoothing
[m,r]=size(vect);
[mbis,d]=size(U);
if m~=mbis
    display('error');
end
kv=zeros(1,m);
for k=1:m
    kv(1,k)=kval(u0-U(k,:),sigma);
end
valf=kv*vect/(kv*ones(m,1));
end

function [val]=kval(u1,sig)
val=exp(-(norm(u1,2)^2)/(2*sig^2));
end