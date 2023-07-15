function FEMsurfaceplot(uh,T,k)

% FEMsurfaceplot(uh,T,k)
% 
% Input:
%       uh   : scalar FEM solution, used for coloring the plot
%       T    : enhanced tetrahedrization 
%       k    : polynomial degree
% 
% Output:
%       The result of this function is a plot on a k-refinement of the
%       surface of the mesh contained in T colored by the function uh
% 
% Last Modified: December 2, 2017

pt = [];
ind =1;
for i = 0:k
    pt = [pt; ...
          (0:k-i)', i*ones(k-i+1,1)];
    list{i+1} = ind+(0:k-i);
    ind = ind+k+1-i; % index number for beginning of next row
end

pt = pt/k;
pt = [1-pt(:,1)-pt(:,2), pt];

ell = k+1:-1:1; % ell(i) = length(list{i})

Tri =[];
for i = 1:k
    t = list{i}(1:ell(i)-1); % bottom left vertex
    Tri = [Tri, [t;1+t;ell(i) + t]];
    t(end) = [];
    Tri = [Tri, [t+1; 1+ell(i)+t; ell(i)+t]];
end

LT = T.faces(1:3, T.faces(4,:)~=0);
X = T.coordinates(1,:); X=pt*X(LT);
Y = T.coordinates(2,:); Y=pt*Y(LT);
Z = T.coordinates(3,:); Z=pt*Z(LT);

d2 = nchoosek(k+2,2);
TRI = bsxfun(@plus, d2*(0:size(LT,2)-1), Tri(:));
TRI = reshape(TRI, 3, k^2*size(LT,2));

faceList = find(T.faces(4,:)~=0);
DOF=computeBDDOF3D(T,k,faceList);

P=bernstein2D(pt(:,2),pt(:,3),k);
U=P*uh(DOF);
trisurf(TRI', X(:), Y(:), Z(:),U(:))
xlabel('x');ylabel('y');
shading interp
