function [P,indices] = bernstein2D(x,y,k)

% [P,ind] = bernstein2D(x,y,k)
%
% Input:
%    [x,y]  : two M x 1 vectors with coordinates in the reference element
%    k      : polynomial degree
% Output:
%    P      : M x dk matrix with Pk basis evaluated at points
%    ind    : 3d array to find where three-index BB basis is located at
% Last modified: April 3, 2015

lambda1=bsxfun(@power,1-x-y,0:k);
lambda2=bsxfun(@power,x,0:k);
lambda3=bsxfun(@power,y,0:k);

indices = location(k);
ind = @(i1,i2,i3) indices(i1+1,i2+1,i3+1);

P=zeros(size(x,1),((k+1)*(k+2))/2);
for a1=0:k
    for a2=0:k-a1
        a3=k-a1-a2;
        P(:,ind(a1,a2,a3))=multinomial(k,a1,a2,a3)*...
                    lambda1(:,a1+1).*lambda2(:,a2+1).*lambda3(:,a3+1);
    end
end

return

function matrix=location(k)
% 3D array with order indices for the BB basis in local FE form

matrix=zeros(k+1,k+1,k+1);
if k==0
    matrix=[1];
    return
end
matrix(k+1,1,1)=1;
matrix(1,k+1,1)=2;
matrix(1,1,k+1)=3;
for j=1:k-1
    matrix(k-j+1,j+1,1)=3+j;
    matrix(1,k-j+1,j+1)=3+k-1+j;
    matrix(j+1,1,k-j+1)=3+2*(k-1)+j;
end
index=3+3*(k-1);
for i=1:k-2
    for j=1:k-1-i
        index=index+1;
        matrix(k-i-j+1,i+1,j+1)=index;
    end
end
return

function y=multinomial(k,i1,i2,i3)
% Bad implementation of multinomial coefficients, onlu valid for small k
y=factorial(k)/(factorial(i1)*factorial(i2)*factorial(i3));
return
