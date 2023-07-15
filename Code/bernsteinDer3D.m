function [Px,Py,Pz]=bernsteinDer3D(x,y,z,k)

% [Px,Py,Pz]=bernsteinDer3D(x,y,z,k)
%
% Input:
%    [x,y,z]  : three M x 1 vectors with coordinates in the reference element
%    k        : polynomial degree
% Output:
%    Px,Py,Pz : M x dk matrices with gradient of Pk basis evaluated at points
% Last modified: April 9, 2015

if k==1
    dim=size(x);
    Px=[-ones(dim) ones(dim) zeros(dim) zeros(dim)];
    Py=[-ones(dim) zeros(dim) ones(dim) zeros(dim)];
    Pz=[-ones(dim) zeros(dim) zeros(dim) ones(dim)];
    return
end

[P,locm1]=bernstein3D(x,y,z,k-1);
[~,loc]  =bernstein3D(0,0,0,k);
indm1 = @(i1,i2,i3,i4) locm1(i1+1,i2+1,i3+1,i4+1);
ind   = @(i1,i2,i3,i4) loc(i1+1,i2+1,i3+1,i4+1);

PxA=zeros(size(x,1),((k+1)*(k+2)*(k+3))/6); 
PxB=PxA; 
PyA=PxA;
PyB=PxA;
PzA=PxA;
PzB=PzA;

for a1=1:k
    for a2=0:k-a1
        for a3=0:k-a1-a2;
            a4=k-a1-a2-a3;
            PxA(:,ind(a1,a2,a3,a4))=-P(:,indm1(a1-1,a2,a3,a4));
            PyA(:,ind(a1,a2,a3,a4))=-P(:,indm1(a1-1,a2,a3,a4));
            PzA(:,ind(a1,a2,a3,a4))=-P(:,indm1(a1-1,a2,a3,a4));
        end
    end
end

for a1=0:k
    for a2=1:k-a1
        for a3=0:k-a1-a2;
            a4=k-a1-a2-a3;
            PxB(:,ind(a1,a2,a3,a4))=P(:,indm1(a1,a2-1,a3,a4));
        end
    end
end

for a1=0:k
    for a2=0:k-a1
        for a3=1:k-a1-a2;
            a4 = k-a1-a2-a3;
                PyB(:,ind(a1,a2,a3,a4))=P(:,indm1(a1,a2,a3-1,a4));
        end
    end
end

for a1=0:k
    for a2=0:k-a1
        for a3=0:k-a1-a2;
            a4 = k-a1-a2-a3;
            if a4>0
                PzB(:,ind(a1,a2,a3,a4))=P(:,indm1(a1,a2,a3,a4-1));
            end
        end
    end
end

Px=k*(PxA+PxB);
Py=k*(PyA+PyB);
Pz=k*(PzA+PzB);
return


