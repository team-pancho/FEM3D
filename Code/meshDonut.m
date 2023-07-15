
function [T,Q] = meshDonut(lev,varargin)

% [T,Q] = meshDonut(lev)
% [T,Q] = meshDonut(lev,Dirichlet)
%
% Input:
%     lev       : refinement level, 1 for no refinement
%     Dirichlet : vector with indices for Dirichlet faces
%                    1 (top), 2 (bottom), 3 (inside), 4 (outside)
% Output:
%     T   : tetrahedron partition of a polyhedral donut
%           default is all faces Dirichlet
%     Q   : quad partition of a polyhedral donut
%           default is all faces Dirichlet
% Last modified: June 28, 2016

% initialize parameter
N= 4; index= 1:N^2; delta= 1;

% generate coordinates
xinitial= -3/2;
yinitial= -3/2;
zinitial= -1/2;
x= (mod(index-1,N))*delta + xinitial;
y= (floor((index-1)/N))*delta + yinitial;
z= ones(1,N^2) * zinitial;
A1= [x' y' z'];
A2= A1; A2(:,3)= delta + A2(:,3);
Q.coordinates= [A1', A2'];

% generate elements
cubes_first_row= [1 2 3 7 11 10 9 5];
cube_matrix_upper_half= [cubes_first_row; cubes_first_row+1; ...
    cubes_first_row+5; cubes_first_row+4];
cube_matrix_lower_half= cube_matrix_upper_half + 16;
Q.elements= [cube_matrix_upper_half;cube_matrix_lower_half];

% generate boundary face
Q.bdface{1}= cube_matrix_lower_half;
Q.bdtag{1} = 'top';

Q.bdface{2}= cube_matrix_upper_half;
Q.bdface{2}(:,:) = Q.bdface{2}([1 4 3 2],:);
Q.bdtag{2} = 'bottom';

externalface_first_row= [1 2 3 4 8 12 16 15 14 13 9 5];
externalface_second_row= [externalface_first_row(2:12) 1];
externalface_matrix= [externalface_first_row; externalface_second_row;...
    externalface_second_row+16; externalface_first_row+16];

Q.bdface{4}=externalface_matrix;
Q.bdface{4}(:,7:12) = Q.bdface{4}([4 1 2 3],7:12);
Q.bdtag{4} = 'outside';

internalface_first_row= [6 7 11 10];
internalface_last_row= [internalface_first_row(2:4) 6];
internalface_matrix= [internalface_first_row; internalface_first_row+16;...
    internalface_last_row+16;internalface_last_row];
Q.bdface{3}= internalface_matrix;
Q.bdface{3}(:,3:4) = Q.bdface{3}([4 1 2 3],3:4);
Q.bdtag{3} = 'inside';

% turn outward to inward
Q.bdface{1}(:,:)= Q.bdface{1}([1 4 3 2],:);
Q.bdface{2}(:,:)= Q.bdface{2}([1 4 3 2],:);
Q.bdface{3}(:,:)= Q.bdface{3}([1 4 3 2],:);
Q.bdface{4}(:,:)= Q.bdface{4}([1 4 3 2],:);

if nargin==2
    listDir=varargin{1};
    listNeu=1:4;
    listNeu(listDir)=[];
else
    listDir=1:4;
    listNeu=[];
end
Q.dirichlet=[];
Q.neumann=[];
for k=listDir
    Q.dirichlet=[Q.dirichlet Q.bdface{k}];
end
for k=listNeu
    Q.neumann=[Q.neumann Q.bdface{k}];
end

T = quad2tet(Q);
for i=1:lev-1
    T = meshTetRefinement(T);
    Q = meshQuadRefinement(Q);
end

end
