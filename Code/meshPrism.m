function T=meshPrism(Nx,Ny,Nz,varargin)

%T=meshPrism(Nx,Ny,Nz)
%T=meshPrism(Nx,Ny,Nz,listDir)
% 
% input
%    Nx  : number of partitions in the x direction
%    Ny  : number of partitions in the y direction
%    Nz  : number of partitions in the z direction
%    listDir : row vector with Dir BC {1:bottom, 2:top, 3:left(x=0), 4
%              4:front(x+y<=1), 5:back(y=0)}.
%              default listDir=[1 2];
% output
%    T  : tetrahedrization of a prism
%
% Last Modified: March 18, 2016

[~,Q] = meshCube(Nx,Ny,Nz);

L = min(Q.coordinates(1,:)/2,Q.coordinates(2,:)/2);
Q.coordinates(1,:)=Q.coordinates(1,:)-L; 
Q.coordinates(2,:)=Q.coordinates(2,:)-L;

split = size(Q.bdface{4},2);
Q.bdface{4}=[Q.bdface{4} Q.bdface{6}];
Q.bdface(6)=[];

Q.bdtag{1}='bottom (z=0)';
Q.bdtag{2}='top (z=1)';
Q.bdtag{3}='left (x=0)';
Q.bdtag{4}='front (x+y<=1)';
Q.bdtag{5}='back (y=0)';
Q.bdtag(6)=[];

% Dirichlet/Neumann division

if nargin==3
    listDir=[1 2];
    listNeu=[3 4 5];
else
    listDir=varargin{1};
    listNeu=1:5;
    listNeu(listDir)=[];
end

T.coordinates=Q.coordinates;

% partition of each prism into 6 tetrahedra

tetra=[1 2 3 5;...
    6 3 5 7;...
    6 5 3 2;...
    1 3 4 5;...
    7 8 4 5;...
    7 4 3 5]';
Ncubes=size(Q.elements,2);
T.elements=Q.elements(tetra(:),:);
T.elements=reshape(T.elements,[4,6*Ncubes]);

% partition of each face into 2 triangles

tri1=[1 4 2;2 4 3]';
tri2=[1 4 3;2 1 3]';

% Separated faces
for i=1:2
    Ni=size(Q.bdface{i},2);
    T.bdface{i}=Q.bdface{i}(tri2(:),:);
    T.bdface{i}=reshape(T.bdface{i},[3,2*Ni]);
end

Ni=size([Q.bdface{4}],2);
T.bdface{4}=[Q.bdface{4}(tri2(:),1:split) Q.bdface{4}(tri1(:),split+1:end)];
T.bdface{4}=reshape(T.bdface{4},[3,2*Ni]);

for i=[3 5]
    Ni=size(Q.bdface{i},2);
    T.bdface{i}=Q.bdface{i}(tri1(:),:);
    T.bdface{i}=reshape(T.bdface{i},[3,2*Ni]);
end

% Dirichlet & Neumann faces
T.dirichlet=[];
T.neumann=[];
for i=listDir
    T.dirichlet=[T.dirichlet T.bdface{i}];
end
for i=listNeu
    T.neumann=[T.neumann T.bdface{i}];
end
    
for i=1:5
    T.bdtag{i}=Q.bdtag{i};
end
return