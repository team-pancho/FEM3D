function formula=quadratureFEM(deg,dim)

% form=quadratureFEM(deg,d)
% Input:
%     deg    : degree of precision of the quad formula
%     d      : dimension
% Output:
%     form   : Stroud Quadrature formula on the d-dimensional simplex
%              (with d=2 or 3) - Nquad x (d+2) matrix
%              with barycentric coordinates and weights
%              (weights are not normalized yet)
% Last modified: June 10, 2016

N=ceil((deg+1)/2);
[t3,w3]=gaussJacobi(N,0);
[t2,w2]=gaussJacobi(N,1);
[t1,w1]=gaussJacobi(N,2);
t3=(t3+1)/2;
t2=(t2+1)/2;
t1=(t1+1)/2;
w3=w3/2;
w2=w2/4;
w1=w1/8;

if dim==3
    table=zeros(N,N,N);
    table(:)=1:N^3;
    x=zeros(N^3,1);
    y=x;z=x;w=x;
    for q1=1:N
        for q2=1:N
            for q3=1:N
                x(table(q1,q2,q3))=t1(q1);
                y(table(q1,q2,q3))=(1-t1(q1))*t2(q2);
                z(table(q1,q2,q3))=(1-t1(q1))*(1-t2(q2))*t3(q3);
                w(table(q1,q2,q3))=w1(q1)*w2(q2)*w3(q3);
            end
        end
    end
    formula=[1-x-y-z x y z 6*w];
elseif dim==2
    table=zeros(N,N);
    table(:)=1:N^2;
    x=zeros(N^2,1); 
    y=x;w=x;
    t1=t2;w1=w2;
    t2=t3;w2=w3;
    for q1=1:N
        for q2=1:N
                x(table(q1,q2))=t1(q1);
                y(table(q1,q2))=(1-t1(q1))*t2(q2);
                w(table(q1,q2))=w1(q1)*w2(q2);
        end
    end
    formula=[1-x-y x y 2*w];
end
return
