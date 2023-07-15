function S=stiffnessMatrices3D(varargin)

% Sh=stiffnessMatrices3DNew(kxx,kxy,kxz,kyy,kyz,kzz,T,k);
% Sh=stiffnessMatrices3DNew(kxx,T,k);
%
% Input:
%    kxx,kxy,kxz,kyy,kyz,kzz : vectorized functions of three variables
%    T  : enhanced triangulation
%    k  : polynomial degree
% Output:
%    Sh : 3 x 3 Cell array dim FE x dim FE stiffness matrices
%         (only upper blocks of cell array are non-empty)
% Last modified: September 30, 2016

% Evaluations of coefficients and basis functions

switch nargin
    case 8
        kxx = varargin{1};
        kxy = varargin{2};
        kxz = varargin{3};
        kyy = varargin{4};
        kyz = varargin{5};
        kzz = varargin{6};
        T = varargin{7};
        k = varargin{8};
    case 3
        kxx = varargin{1};
        kxy = kxx; kxz = kxx; kyy = kxx; kyz = kxx; kzz = kxx;
        T = varargin{2};
        k = varargin{3};
end

formula=quadratureFEM(3*k+1,3);
x=T.coordinates(1,:); x=formula(:,1:4)*x(T.elements);
y=T.coordinates(2,:); y=formula(:,1:4)*y(T.elements);
z=T.coordinates(3,:); z=formula(:,1:4)*z(T.elements);

K{1,1}=kxx(x,y,z);
K{1,2}=kxy(x,y,z);
K{1,3}=kxz(x,y,z);
K{2,2}=kyy(x,y,z);
K{2,3}=kyz(x,y,z);
K{3,3}=kzz(x,y,z);
[Px,Py,Pz]=bernsteinDer3D(formula(:,2),formula(:,3),formula(:,4),k);

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

% Loop over quadrature points

dk=size(Px,2); 
Nelt=size(T.elements,2);
for i=1:3
    for j=i:3
        S{i,j}=zeros(dk,dk*Nelt);
    end
end

for q=1:size(formula,1)
    P{1,1}=formula(q,5)*Px(q,:)'*Px(q,:);
    P{1,2}=formula(q,5)*Px(q,:)'*Py(q,:);
    P{1,3}=formula(q,5)*Px(q,:)'*Pz(q,:);
    P{2,1}=P{1,2}';%formula(q,5)*Py(q,:)'*Px(q,:);
    P{2,2}=formula(q,5)*Py(q,:)'*Py(q,:);
    P{2,3}=formula(q,5)*Py(q,:)'*Pz(q,:);
    P{3,1}=P{1,3}';%formula(q,5)*Pz(q,:)'*Px(q,:);
    P{3,2}=P{2,3}';%formula(q,5)*Pz(q,:)'*Py(q,:);
    P{3,3}=formula(q,5)*Pz(q,:)'*Pz(q,:);
    
    for i=1:3
        for j=i:3
            for l=1:3
                for m=1:3
                    S{i,j}=S{i,j}+kron(c{l,j}.*c{m,i}.*K{i,j}(q,:),P{m,l});
                end
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