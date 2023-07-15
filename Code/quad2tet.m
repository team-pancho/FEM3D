function T=quad2tet(Q)

% T=quad2tet(Q)
% 
% Input:
%    Q : basic quad partition
% Output:
%    T : tetrahedral partition with the same nodes and boundary 
% Last modified: July 15, 2015

T.coordinates=Q.coordinates;

% partition of each quad into 6 tetrahedra

tetra=[1 2 4 5;...
    6 4 5 8;...
    6 5 4 2;...
    2 3 4 6;...
    8 7 6 3;...
    3 8 4 6]';
Ncubes=size(Q.elements,2);
T.elements=Q.elements(tetra(:),:);
T.elements=reshape(T.elements,[4,6*Ncubes]);

% partition of each face into 2 triangles

tri=[1 4 2;2 4 3]';

% Dirichlet faces

NDir=size(Q.dirichlet,2);
if isempty(Q.dirichlet)
    T.dirichlet=[];
else
    T.dirichlet=Q.dirichlet(tri(:),:);
    T.dirichlet=reshape(T.dirichlet,[3,2*NDir]);
end
% Neumann faces

NNeu=size(Q.neumann,2);
if isempty(Q.neumann)
    T.neumann=[];
else
    T.neumann=Q.neumann(tri(:),:);
    T.neumann=reshape(T.neumann,[3,2*NNeu]);
end
    
% Separated face

if isfield(Q,'bdface')
    for i=1:length(Q.bdface)
        Ni=size(Q.bdface{i},2);
        T.bdface{i}=Q.bdface{i}(tri(:),:);
        T.bdface{i}=reshape(T.bdface{i},[3,2*Ni]);
    end
    T.bdtag=Q.bdtag;
end

return
