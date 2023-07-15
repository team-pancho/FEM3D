function S=stiffnessMatricesCC3D(T,k)

% Sh=stiffnessMatricesCC3D(T,k);
%
% Input:
%    T  : enhanced triangulation
%    k  : polynomial degree
% Output:
%    Sh : 3 x 3 Cell array dim FE x dim FE stiffness matrices
%         (only upper blocks of cell array are non-empty)
% Last modified: September 30, 2016

% Matrices on the reference element

formula=quadratureFEM(2*k-1,3);

[Px,Py,Pz]=bernsteinDer3D(formula(:,2),formula(:,3),formula(:,4),k);
wPx=bsxfun(@times,formula(:,5),Px);
wPy=bsxfun(@times,formula(:,5),Py);
wPz=bsxfun(@times,formula(:,5),Pz);
P{1,1}=Px'*wPx;
P{1,2}=Px'*wPy;
P{1,3}=Px'*wPz;
P{2,1}=P{1,2}';
P{2,2}=Py'*wPy;
P{2,3}=Py'*wPz;
P{3,1}=P{1,3}';
P{3,2}=P{2,3}';
P{3,3}=Pz'*wPz;

% Geometric coefficients for change of variables

x12=T.coordinates(1,T.elements(2,:))-T.coordinates(1,T.elements(1,:));  
y12=T.coordinates(2,T.elements(2,:))-T.coordinates(2,T.elements(1,:));  
z12=T.coordinates(3,T.elements(2,:))-T.coordinates(3,T.elements(1,:)); 
x13=T.coordinates(1,T.elements(3,:))-T.coordinates(1,T.elements(1,:));   
y13=T.coordinates(2,T.elements(3,:))-T.coordinates(2,T.elements(1,:));
z13=T.coordinates(3,T.elements(3,:))-T.coordinates(3,T.elements(1,:));
x14=T.coordinates(1,T.elements(4,:))-T.coordinates(1,T.elements(1,:));   
y14=T.coordinates(2,T.elements(4,:))-T.coordinates(2,T.elements(1,:));
z14=T.coordinates(3,T.elements(4,:))-T.coordinates(3,T.elements(1,:));

% Entries of CK

sqdet=1./(6*sqrt(T.volume));
c{1,1}=sqdet.*(y13.*z14 - y14.*z13);
c{1,2}=sqdet.*(x14.*z13 - x13.*z14);
c{1,3}=sqdet.*(x13.*y14 - x14.*y13);
c{2,1}=sqdet.*(y14.*z12 - y12.*z14);
c{2,2}=sqdet.*(x12.*z14 - x14.*z12);
c{2,3}=sqdet.*(x14.*y12 - x12.*y14);
c{3,1}=sqdet.*(y12.*z13 - z12.*y13);
c{3,2}=sqdet.*(x13.*z12 - x12.*z13);
c{3,3}=sqdet.*(x12.*y13 - x13.*y12);

% Matrices on the physical elements

dk=size(Px,2); 
Nelt=size(T.elements,2);
for i=1:3
    for j=i:3
        S{i,j}=zeros(dk,dk*Nelt);
        for el=1:3
            for m=1:3
                S{i,j}=S{i,j}+kron(c{el,j}.*c{m,i},P{m,el});
            end
        end
    end
end

% Assembly 

[~,Cols]=DOF3D(T,k);
Rows=permute(Cols,[2 1 3]);
for i=1:3
    for j=i:3
        S{i,j}=sparse(Rows(:),Cols(:),S{i,j}(:));
    end
end
return