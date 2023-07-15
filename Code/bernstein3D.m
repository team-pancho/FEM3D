    function [P,indices]=bernstein3D(x,y,z,k)

% [P,ind]=bernstein3D(x,y,z,k)
%
% Input:
%    [x,y,z]: three M x 1 vectors with coordinates in the reference element
%    k      : polynomial degree
% Output:
%    P      : M x dk matrix with Pk basis evaluated at points
%    ind    : 4d array to find where four-index BB basis is located at
% Last modified: April 9, 2015

if k==0
    P=ones(size(x));
    indices=1;
    return
end

lambda1=bsxfun(@power,1-x-y-z,0:k);
lambda2=bsxfun(@power,x,0:k);
lambda3=bsxfun(@power,y,0:k);
lambda4=bsxfun(@power,z,0:k);

indices = location(k);
ind = @(i1,i2,i3,i4) indices(i1+1,i2+1,i3+1,i4+1);

P=zeros(size(x,1),((k+1)*(k+2)*(k+3))/6);
for a1=0:k
    for a2=0:k-a1
        for a3=0:k-a1-a2
            a4=k-a1-a2-a3;
            P(:,ind(a1,a2,a3,a4))=multinomial(k,a1,a2,a3,a4)*...
                    lambda1(:,a1+1).*lambda2(:,a2+1).*lambda3(:,a3+1).*lambda4(:,a4+1);
        end
    end
end

return

function matrix=location(k)
% 4D array with order indices for the BB basis in local FE form

matrix=zeros(k+1,k+1,k+1,k+1);

matrix(k+1,1,1,1)=1;
matrix(1,k+1,1,1)=2;
matrix(1,1,k+1,1)=3;
matrix(1,1,1,k+1)=4;

for j=1:k-1
    matrix(k-j+1,j+1,1,1)=4+j;
    matrix(1,k-j+1,j+1,1)=4+k-1+j;
    matrix(k-j+1,1,j+1,1)=4+2*(k-1)+j;
    matrix(1,k-j+1,1,j+1)=4+3*(k-1)+j;
    matrix(1,1,k-j+1,j+1)=4+4*(k-1)+j;
    matrix(k-j+1,1,1,j+1)=4+5*(k-1)+j;
end

sofar=4+6*(k-1);
dimfac=(k-1)*(k-2)/2;

for i=1:k-2
    di=(i-1)*(k-2)-(i-2)*(i-1)/2;
    for j=1:k-1-i
        matrix(k-i-j+1,i+1,j+1,1)=sofar+di+j;
        matrix(k-i-j+1,i+1,1,j+1)=sofar+dimfac+di+j;
        matrix(k-i-j+1,1,i+1,j+1)=sofar+2*dimfac+di+j;
        matrix(1,i+1,j+1,k-i-j+1)=sofar+3*dimfac+di+j;
    end
end

index = 4+6*(k-1)+4*dimfac;
for p=1:k-3
    for i=1:k-2-p
        for j=1:k-1-p-i
            index=index+1;
            matrix(k-i-j-p+1,p+1,i+1,j+1) = index;
        end
    end
end

return


function y=multinomial(k,i1,i2,i3,i4)
% Bad implementation of multinomial coefficients, only valid for small k
y=factorial(k)/(factorial(i1)*factorial(i2)*factorial(i3)*factorial(i4));
return