function Ch=convectionMatrices3D(bx, by, bz,T,k)

% Ch=convectionMatrices3D(bx, by, bz,T,k);
%
% Input:
%    bx,by,bz : vectorized functions of three variables
%    T  : enhanced triangulation
%    k  : polynomial degree
% Output:
%    Ch : 1 x 3 Cell array dim FE x dim FE convection matrices
%         \int b_i P_a \partial_i P_b
%           
% Last modified: June 3, 2016

% Evaluations of coefficients and basis functions

formula=quadratureFEM(3*k+1,3);
x=T.coordinates(1,:); x=formula(:,1:4)*x(T.elements);
y=T.coordinates(2,:); y=formula(:,1:4)*y(T.elements);
z=T.coordinates(3,:); z=formula(:,1:4)*z(T.elements);

B{1}=bx(x,y,z);
B{2}=by(x,y,z);
B{3}=bz(x,y,z);
P=bernstein3D(formula(:,2),formula(:,3),formula(:,4),k);
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

d{1,1}=1/6.*(y13.*z14 - y14.*z13);
d{1,2}=1/6.*(x14.*z13 - x13.*z14);
d{1,3}=1/6.*(x13.*y14 - x14.*y13);
d{2,1}=1/6.*(y14.*z12 - y12.*z14);
d{2,2}=1/6.*(x12.*z14 - x14.*z12);
d{2,3}=1/6.*(x14.*y12 - x12.*y14);
d{3,1}=1/6.*(y12.*z13 - z12.*y13);
d{3,2}=1/6.*(x13.*z12 - x12.*z13);
d{3,3}=1/6.*(x12.*y13 - x13.*y12);

% Loop over quadrature points

dk=size(Px,2); 
Nelt=size(T.elements,2);
for i=1:3
    Ch{i}=zeros(dk,dk*Nelt);
end

for q=1:size(formula,1)
    Q{1}=formula(q,5)*P(q,:)'*Px(q,:);
    Q{2}=formula(q,5)*P(q,:)'*Py(q,:);
    Q{3}=formula(q,5)*P(q,:)'*Pz(q,:);
    
    for i=1:3
        for j=1:3
            Ch{i}=Ch{i}+kron(d{j,i}.*B{i}(q,:),Q{j});
        end
    end
end


% Assembly 
[~,Cols]=DOF3D(T,k);
Rows=permute(Cols,[2 1 3]);
for i=1:3
    Ch{i}=sparse(Rows(:),Cols(:),Ch{i}(:));
end
return