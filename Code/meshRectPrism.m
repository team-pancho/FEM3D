function [T,Q]=meshRectPrism(dx,dy,dz,lev,varargin)

% [T,Q]=meshRectPrism(dx,dy,dz,lev)
% [T,Q]=meshRectPrism(dx,dy,dz,lev,listDir)
% Input
%    dx  : dimension of the prism in the x direction
%    dy  : dimension of the prism in the y direction
%    dz  : dimension of the prism in the z direction
%    lev : nof refinements per unit dimension
%    listDir : list of Dirichlet faces: bottom(1), top(2), left(3),
%              right(4), back(5), front(6)
% Output
%    T   : tetrahedral partition data structure
%          default: the top and bottom faces are Dirichlet
%    Q   : quadrilateral mesh of the same data structure
%          default: the top and bottom faces are Diriclet
% Last modified: January 10, 2019

Nx=lev*dx;Ny=lev*dy;Nz=lev*dz;
Nx=Nx+1;
Ny=Ny+1;
Nz=Nz+1;

% creating the cubic partition

list=reshape(1:Nx*Ny*Nz,Nx,Ny,Nz);
c=list(1:Nx-1,1:Ny-1,1:Nz-1);
q=[c(:) c(:)+1 c(:)+Nx+1 c(:)+Nx];
Q.elements=[q q+Nx*Ny]';

% coordinates of the nodes
% you can change the cube to another paralellepiped here

x=linspace(0,dx,Nx);        
y=linspace(0,dy,Ny);
z=linspace(0,dz,Nz);
[y,x,z]=meshgrid(y,x,z);
Q.coordinates=[x(:) y(:) z(:)]';

reverse=[1 4 3 2];

% faces 1 (bottom z=0) & 2 (top z=1)
c=list(1:Nx-1,1:Ny-1,1);
face{1}=[c(:) c(:)+1 c(:)+Nx+1 c(:)+Nx];
face{2}=(Nz-1)*Nx*Ny+face{1};
face{2}=face{2}(:,reverse);

% faces 3 (left y=0) & 4 (right y=1)
c=list(1:Nx-1,1,1:Nz-1);
face{3}=[c(:) c(:)+Nx*Ny c(:)+Nx*Ny+1 c(:)+1];
face{4}=(Ny-1)*Nx+face{3};
face{4}=face{4}(:,reverse);

% faces 5 (back x=0) & 6 (front x=1)
c=list(1,1:Ny-1,1:Nz-1);
face{5}=[c(:) c(:)+Nx c(:)+Nx*Ny+Nx c(:)+Nx*Ny];
face{6}=(Nx-1)+face{5};
face{6}=face{6}(:,reverse);

for i=1:6
    Q.bdface{i}=face{i}';
end
Q.bdtag{1}='bottom (z=0)';
Q.bdtag{2}=['top (z=' num2str(Nz-1) ')'];
Q.bdtag{3}='left (y=0)';
Q.bdtag{4}=['right (y=' num2str(Ny-1) ')'];
Q.bdtag{5}='back (x=0)';
Q.bdtag{6}=['front (x=' num2str(Nx-1) ')'];

% Dirichlet/Neumann division

if nargin==4
    listDir=[1 2];
    listNeu=[3 4 5 6];
else
    listDir=varargin{1};
    listNeu=1:6;
    listNeu(listDir)=[];
end
Q.dirichlet=[];
Q.neumann=[];
for i=listDir
    Q.dirichlet=[Q.dirichlet face{i}'];
end
for i=listNeu
    Q.neumann=[Q.neumann face{i}'];
end
T = quad2tet(Q);
return    

