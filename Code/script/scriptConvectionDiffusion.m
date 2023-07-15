% BENCHMARK SCRIPT
% This script tests the equation -div(k grad u)+ b grad u + c u = f
% with mixed boundary conditions
% It should always be tested with pure Dirichlet and pure Neumann BC to be
% sure that everything works
% Last Modified: 9/30/16

clear; clc; close all;

u     = @(x,y,z) cos(pi*x).*sin(pi*y).*cos(pi*z);
ux    = @(x,y,z) -pi*sin(pi*x).*sin(pi*y).*cos(pi*z);
uy    = @(x,y,z) pi*cos(pi*x).*cos(pi*y).*cos(pi*z);
uz    = @(x,y,z) -pi*cos(pi*x).*sin(pi*y).*sin(pi*z);
uxx   = @(x,y,z) -pi^2*u(x,y,z);
uxy   = @(x,y,z) -pi^2*sin(pi*x).*cos(pi*y).*cos(pi*z);
uxz   = @(x,y,z) pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z);
uyy   = @(x,y,z) -pi^2*u(x,y,z);
uyz   = @(x,y,z) -pi^2*cos(pi*x).*cos(pi*y).*sin(pi*z);
uzz   = @(x,y,z) -pi^2*u(x,y,z);

c=@(x,y,z) 1+x.^2+y.^2+z.^2;

% Constant Kappa

% k11=@(x,y,z) 3+0*x; 
% k22=@(x,y,z) 5+0*x; 
% k33=@(x,y,z) 6+0*x; 
% k12=@(x,y,z) 0 +0*x; 
% k13=@(x,y,z) 0 +0*x;
% k23=@(x,y,z) 0 +0*x;
% k21=@(x,y,z) k12(x,y,z);
% k31=@(x,y,z) k13(x,y,z);
% k32=@(x,y,z) k23(x,y,z);
% 
% k11x = @(x,y,z) 0+0*x;
% k22y = @(x,y,z) 0+0*x;
% k33z = @(x,y,z) 0+0*x;
% k12x = @(x,y,z) 0+0*x;
% k12y = @(x,y,z) 0+0*x;
% k12z = @(x,y,z) 0+0*x;
% k13x = @(x,y,z) k12x(x,y,z);
% k23y = @(x,y,z) k12y(x,y,z);
% k21y = @(x,y,z) k12y(x,y,z);
% k31z = @(x,y,z) k12z(x,y,z);
% k32z = @(x,y,z) k12z(x,y,z);

% Variable Kappa

k11=@(x,y,z) 2./(1+x.^2+y.^2+z.^2)+1+x.^2;
k22=@(x,y,z) 2./(1+x.^2+y.^2+z.^2)+1+y.^2;
k33=@(x,y,z) 2./(1+x.^2+y.^2+z.^2)+1+z.^2;
k12=@(x,y,z) 1./(1+x.^2+y.^2+z.^2);
k13 = @(x,y,z) k12(x,y,z);
k23 = @(x,y,z) k12(x,y,z);

k11x = @(x,y,z) -4*x./(1+x.^2+y.^2+z.^2).^2+2*x;
k22y = @(x,y,z) -4*y./(1+x.^2+y.^2+z.^2).^2+2*y;
k33z = @(x,y,z) -4*z./(1+x.^2+y.^2+z.^2).^2+2*z;
k12x = @(x,y,z) -2*x./(1+x.^2+y.^2+z.^2).^2;
k12y = @(x,y,z) -2*y./(1+x.^2+y.^2+z.^2).^2;
k12z = @(x,y,z) -2*z./(1+x.^2+y.^2+z.^2).^2;
k13x = @(x,y,z) k12x(x,y,z);
k23y = @(x,y,z) k12y(x,y,z);
k13z = @(x,y,z) k12z(x,y,z);
k23z = @(x,y,z) k12z(x,y,z);

bx = @(x,y,z) sin(y.*z);
by = @(x,y,z) z.*exp(x.*z);
bz = @(x,y,z) 3.*(x.*y).^2;

g1=@(x,y,z) k11(x,y,z).*ux(x,y,z)+k12(x,y,z).*uy(x,y,z)+k13(x,y,z).*uz(x,y,z);
g2=@(x,y,z) k12(x,y,z).*ux(x,y,z)+k22(x,y,z).*uy(x,y,z)+k23(x,y,z).*uz(x,y,z);
g3=@(x,y,z) k13(x,y,z).*ux(x,y,z)+k23(x,y,z).*uy(x,y,z)+k33(x,y,z).*uz(x,y,z);

g={g1,g2,g3};

f = @(x,y,z) c(x,y,z).*u(x,y,z) ...
             + bx(x,y,z).*ux(x,y,z)...
             + by(x,y,z).*uy(x,y,z)...
             + bz(x,y,z).*uz(x,y,z)...
             - k11x(x,y,z).*ux(x,y,z) - k11(x,y,z).*uxx(x,y,z)...
             - k12x(x,y,z).*uy(x,y,z) - k12(x,y,z).*uxy(x,y,z)...
             - k13x(x,y,z).*uz(x,y,z) - k13(x,y,z).*uxz(x,y,z)...
             - k12y(x,y,z).*ux(x,y,z) - k12(x,y,z).*uxy(x,y,z)...
             - k22y(x,y,z).*uy(x,y,z) - k22(x,y,z).*uyy(x,y,z)...
             - k23y(x,y,z).*uz(x,y,z) - k23(x,y,z).*uyz(x,y,z)...
             - k13z(x,y,z).*ux(x,y,z) - k13(x,y,z).*uxz(x,y,z)...
             - k23z(x,y,z).*uy(x,y,z) - k23(x,y,z).*uyz(x,y,z)...
             - k33z(x,y,z).*uz(x,y,z) - k33(x,y,z).*uzz(x,y,z);
         
% discrete parameters
k=input('Polynomial degree: ');
BCs=input('Dirichlet (1) - Neumann (2) - Mixed (3): ');
switch BCs
    case 1
        DBC=1:6;
    case 2
        DBC=[];
    case 3
        DBC=[1 2 4];
end
for j=1:6
    % mesh generation and post-processing
    T=meshCube(j,j,j,[DBC]); 
    T=edgesAndFaces(T); 
    T=enhanceGrid3D(T);
    % build matrices and rhs vectors
    S=stiffnessMatrices3D(k11,k12,k13,k22,k23,k33,T,k);
    S=S{1,1}+S{2,2}+S{3,3}+S{1,2}+S{1,2}'+S{1,3}+S{1,3}'+S{2,3}+S{2,3}';
    M=massMatrix3D(c,T,k);
    C = convectionMatrices3D(bx,by,bz,T,k);
    C = C{1} + C{2} + C{3};
    RHS = loadVector3D(f,T,k)+neumannBC3D(g,T,k);
    [uh,Dir,Free] = dirichletBC3D(u,T,k);
    % FEM Solution
    RHS = RHS - (S(:,Dir) + M(:,Dir) + C(:,Dir))*uh(Dir);
    uh(Free) =(S(Free,Free)+M(Free,Free) + C(Free,Free))\RHS(Free);
    % computation of errors
    [eh(j),hh(j)]=errorFEM3D({u,ux,uy,uz},uh,T,k);
    h(j) = 1/j;
end

% H1 errors
loglog(h,hh, '-or');
hold on
loglog(h,h.^k,'-r');

% L2 errors
loglog(h, eh, '-ob');
hold on
loglog(h,h.^(k+1),'-b');