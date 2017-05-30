function [D1,D2]=grad2D(f,T,U1,U2,step)
EPS=0.0001;
fact=10/step;
l=size(T,1);
[m]=size(U1,1);
D1=zeros(m,1);
D2=zeros(m,1);
for i1=2:l
    u1=T(i1);
    for i2=2:l
        u2=T(i2);
        j=find((fact*U1 + U2)==(fact*u1+u2));
        jprecx=find(abs((fact*U1 + U2)-(fact*(u1-step)+u2))<EPS);
        jprecy=find(abs((fact*U1 + U2)-(fact*u1+(u2-step)))<EPS);
        if (length(jprecx)~=1) | (length(jprecy)~=1)
           error('Problem in grad2D');
        end
        D1(j)=(f(j)-f(jprecx))/step;
        D2(j)=(f(j)-f(jprecy))/step;
    end
end
end