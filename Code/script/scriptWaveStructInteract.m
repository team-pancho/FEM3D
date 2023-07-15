% Script to solve the wave-structure interaction problem
% using FEM in both interior and exterior domain with an 
% impedence boundary condition on the exterior boundary
% 
% Last Modified: January 24, 2017


clear;  close all;

% physical and wave parameters
rho = 3.4;
lam = 3;
mu = 2;
e = [1,1,1];
d = 1/sqrt(3)*e;
option = 'p'; %change this based on your definition of e and d
ad = [1,0,0];

% material and input functions
f = @(x,y,z,s) 0*x;
rhof  = @(x,y,z) 0*x + rho;
lamf  = @(x,y,z) 0*x + lam;
muf  = @(x,y,z) 0*x + mu;

% if e*d.' == 0 
%     c = sqrt()

% frequency
s = 2;

[U,V,W,uInc,g,G] = solutionWaveStructure(d,e,rho,lam,mu,ad,s,option);

% mesh and poly deg.
k = input('polynomial degree: ');    

xmin = 0.33;   % approximate boundary of inner cube
xmax = 0.66;
ymin = 0.33;
ymax = 0.66;
zmin = 0.33;
zmax = 0.66;

planes = [ 1  0  0 -xmin;...
          -1  0  0  xmax;...
           0  1  0 -ymin;...
           0 -1  0  ymax;...
           0  0  1 -zmin;...
           0  0 -1  zmax];

% Cube with the interior removed
n = input('refinement level: ');
TT = meshCube(n,n,n,[]);
TT = edgesAndFaces(TT);
TT = enhanceGrid3D(TT);
eltsm = planeCut(TT,planes);
totalElts = 1:size(TT.elements,2);
eltsp = setdiff(totalElts,eltsm);
[Tm, subTm] = subtet(TT,eltsm);
[Tp, subTp] = subtet(TT,eltsp);
em = embed(TT,subTm,k);
ep = embed(TT,subTp,k);


[Mh,b] = couplingMatrices(TT,Tp,Tm,ep,em,rhof,lamf,muf,k,s,g,G,uInc);
uh = Mh\b;

dim = dimFEMspace(Tm,k);
% dimExt = dimFEMspace(Tp,k);
vh = uh(dim+1:2*dim);
wh = uh(2*dim+1:3*dim);
awh = uh(3*dim+1:end);
uh(dim+1:end) = [];

[ehu,hhu]=errorFEM3D(U,uh,Tm,k);
[ehv,hhv]=errorFEM3D(V,vh,Tm,k);
[ehw,hhw]=errorFEM3D(W,wh,Tm,k);
[ehaw,hhaw]=errorFEM3D(uInc,awh,Tp,k);

ErrElL2 = ehu+ehv+ehw;
ErrElH1 = hhu+hhv+hhw;
ErrL2=ehaw;
ErrH1=hhaw;

disp(ErrElL2)
disp(ErrElH1)
disp(ErrL2)
disp(ErrH1)