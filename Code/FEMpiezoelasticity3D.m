function [u1h,u2h,u3h,psi] = FEMpiezoelasticity3D(mu,lam,rho,e,kap,f,uD,...
                                                  sigmaeta,T,k,varargin)
%
% [u1h,u2h,u3h,psi] = FEMpiezoElasticity3D(mu,lam,rho,e,kap,...
%                      {fx,fy,fz,fpsi},{ux,uy,uz,psiD},...
%                      {sigxx,sigxy,sigxz,sigyy,sigyz,sigzz,etax,etay,...
%                       etaz},T,k)
% [u1h,u2h,u3h,psi] = FEMpiezoElasticity3D(mu,lam,rho,e,kap,...
%                      {fx,fy,fz,fpsi},{ux,uy,uz,psiD},...
%                      {sigxx,sigxy,sigxz,sigyy,sigyz,sigzz,etax,etay,...
%                       etaz},T,k,s)
% [u1h,u2h,u3h,psi] = FEMpiezoElasticity3D(mu,lam,rho,e,kap,...
%                      {fx,fy,fz,fpsi},{ux,uy,uz,psiD},...
%                      {sigxx,sigxy,sigxz,sigyy,sigyz,sigzz,etax,etay,...
%                       etaz},T,k,s,rhsVec)
%
% Input:
% mu,lam           : Vectorized function handles of three variables
% rho              : vectorized function of three variables (density)
% e                : 3 x 6 cell array of vectorized function handles 
%                    for piezoelectric tensor
% kap              : 1 X 6 cell array of vectorized function handles
%                    for dielectric tensor
% fx, fy, fz, fpsi : Vectorized functions of three variables (load vector)
% ux, uy, uz, psiD : Vectorized functions of three variables (Dirichlet BC)
% sigxx,...,etax,..: Nine vectorized functions of three variables, six for
%                    the stress and three for the function used for 
%                    Neumann BC in the last equation 
% T                : Data structure. Basic FEM triangulation
% k                : scalar. Polynomial degree
% s                : Complex wave number (s=-1i k)
% rhsVec           : 3*dimVh-vector of tests for RHS
%
% Output:
% u1h,u2h,u3h      : Ndof x 1 Vectors. Vh coefficients of the displacement
%                    components.
% psi              : Ndof x 1 vector of piezoelectric coefficients
%
% Last modified July 13, 2016.

T = edgesAndFaces(T);
T = enhanceGrid3D(T);

fh =loadVector3D(f,T,k); 
fx = fh(:,1);
fy = fh(:,2);
fz = fh(:,3);
fpsi = fh(:,4);

traction=neumannBC3D({sigmaeta{1:6}},T,k);
tx=traction(:,1);
ty=traction(:,2);
tz=traction(:,3);

ned = neumannBC3D({sigmaeta{7:9}},T,k);

if nargin==12
    rhsVec = varargin{2};
else
    rhsVec = zeros(size([fx;fy;fz;fpsi]));
end

Sm = stiffnessMatrices3D(mu,mu,mu,mu,mu,mu,T,k);
Sl = stiffnessMatrices3D(lam,lam,lam,lam,lam,lam,T,k);

E14 = stiffnessMatricesNonSymmetric3D(e{1,1},e{2,1},e{3,1},e{1,2},...
                               e{2,2},e{3,2},e{1,3},e{2,3},e{3,3},T,k);
E14 = E14{1,1} + E14{1,2} + E14{1,3} + E14{2,1} + E14{2,2} + E14{2,3}...
    + E14{3,1} + E14{3,2} + E14{3,3};
E24 = stiffnessMatricesNonSymmetric3D(e{1,2},e{2,2},e{3,2},e{1,4},...
                               e{2,4},e{3,4},e{1,5},e{2,5},e{3,5},T,k);
E24 = E24{1,1} + E24{1,2} + E24{1,3} + E24{2,1} + E24{2,2} + E24{2,3}...
    + E24{3,1} + E24{3,2} + E24{3,3};
E34 = stiffnessMatricesNonSymmetric3D(e{1,3},e{2,3},e{3,3},e{1,5},...
                               e{2,5},e{3,5},e{1,6},e{2,6},e{3,6},T,k);
E34 = E34{1,1} + E34{1,2} + E34{1,3} + E34{2,1} + E34{2,2} + E34{2,3}...
    + E34{3,1} + E34{3,2} + E34{3,3};

Kap = stiffnessMatrices3D(kap{1},kap{2},kap{3},kap{4},kap{5},kap{6}, T,k);
Kap = Kap{1,1} + Kap{1,2} + Kap{1,2}' + Kap{1,3} + Kap{1,3}' + Kap{2,2}...
    + Kap{2,3} + Kap{2,3}' + Kap{3,3};

% MAT = [(C \varepsilon(w), \varepsilon(w) (e \grad \phi, \grad \phi);...
%       -(e^T \varepsilon(w), \grad \phi)  (\kappa \grad \phi, \grad \phi)]  

MAT=[2*Sm{1,1}+Sl{1,1}+Sm{2,2}+Sm{3,3},Sl{1,2}+Sm{1,2}',...
     Sl{1,3}+Sm{1,3}',                 E14;...
     Sm{1,2}+Sl{1,2}',                 2*Sm{2,2}+Sm{1,1}+Sl{2,2}+Sm{3,3}...
     Sl{2,3}+Sm{2,3}',                 E24;...
     Sl{1,3}'+Sm{1,3},                 Sl{2,3}'+Sm{2,3},...
     2*Sm{3,3}+Sm{1,1}+Sm{2,2}+Sl{3,3},E34;...
     -E14',                            -E24',...
     -E34',                            Kap];
  
if nargin > 10
     s = varargin{1};
     Mh = massMatrix3D(rho,T,k);
     O=sparse(size(Mh,1),size(Mh,2));
     MAT=MAT+s^2*[Mh O O O;...
                  O Mh O O;...
                  O O Mh O;...
                  O O O  O];    
end
 
[uh,~,free] = dirichletBC3D(uD,T,k);
uh=uh(:);
RHS = [fx+tx; fy+ty;fz+tz; fpsi - ned]+rhsVec-MAT*uh;

Ndof = size(fx,1);
free = [free(:); Ndof+free(:);2*Ndof+free(:);3*Ndof + free(:)];
uh(free) = MAT(free,free)\RHS(free);
u1h = reshape(uh,[Ndof,4]);
u2h = u1h(:,2);
u3h = u1h(:,3);
psi = u1h(:,4);
u1h(:,2:4) = [];

end