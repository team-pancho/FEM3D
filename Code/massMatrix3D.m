function Mh=massMatrix3D(c,T,k)

% Mh=massMatrix3D(c,T,k)
%
% Input:
%    c  : vectorized function of three variables
%    T  : enhanced triangulation
%    k  : polynomial degree
% Output:
%    Mh : dim P_k(T_h) x dim P_k(T_h) mass matrix
% Last modified: April 17, 2015

formula=quadratureFEM(4*k+3,3);
x=T.coordinates(1,:); x=formula(:,1:4)*x(T.elements);
y=T.coordinates(2,:); y=formula(:,1:4)*y(T.elements);
z=T.coordinates(3,:); z=formula(:,1:4)*z(T.elements);
c=bsxfun(@times,T.volume,c(x,y,z));
P=bernstein3D(formula(:,2),formula(:,3),formula(:,4),k);

dk=size(P,2); 
Nelt=size(T.elements,2);
Mh=zeros(dk,dk*Nelt);

for q=1:size(formula,1)
    Mh=Mh+kron(c(q,:),formula(q,5)*P(q,:)'*P(q,:));
end
[~,Cols]=DOF3D(T,k);
Rows=permute(Cols,[2 1 3]);
Mh=sparse(Rows(:),Cols(:),Mh(:));

return