function [u1h,u2h,u3h] = FEMelasticity3D(mu,lam,f,uD,sigma,T,k,varargin)
%
% [u1h,u2h,u3h] = FEMelasticity3D(mu,lam,{fx,fy,fz},{ux,uy,uz},...
%                      {sigxx,sigxy,sigxz,sigyy,sigyz,sigzz},T,k)
% [u1h,u2h,u3h] = FEMelasticity3D(mu,lam,{fx,fy,fz},{ux,uy,uz},...
%                      {sigxx,sigxy,sigxz,sigyy,sigyz,sigzz},T,k,rho,s)
% [u1h,u2h,u3h] = FEMelasticity3D(mu,lam,{fx,fy,fz},{ux,uy,uz},...
%                      {sigxx,sigxy,sigxz,...},T,k,rho,s,rhsVec)
%
% Input:
% mu,lam        : Vectorized function handles of three variables
% fx, fy, fz    : Vectorized functions of three variables (load vector)
% ux, uy, uz    : Vectorized functions of three variables (Dirichlet BC)
% sigxx,...     : Six vectorized functions of three variables
% T             : Data structure. Basic FEM triangulation
% k             : scalar. Polynomial degree
% rho           : vectorized function of three variables (density)
% s             : Complex wave number (s=-1i k)
% rhsVec        : 3*dimVh-vector of tests for RHS
%
% Output:
% u1h,u2h,u3h   : Ndof x 1 Vectors. Vh coefficients of the displacement
%                 components.
%
% Last modified December 9, 2015.

T = edgesAndFaces(T);
T = enhanceGrid3D(T);

fh =loadVector3D(f,T,k); 
fx = fh(:,1);
fy = fh(:,2);
fz = fh(:,3);

traction=neumannBC3D(sigma,T,k);
tx=traction(:,1);
ty=traction(:,2);
tz=traction(:,3);

if nargin==10
    rhsVec = varargin{3};
else
    rhsVec = zeros(size([fx;fy;fz]));
end

Sm = stiffnessMatrices3D(mu,mu,mu,mu,mu,mu,T,k);
Sl = stiffnessMatrices3D(lam,lam,lam,lam,lam,lam,T,k);

MAT=[2*Sm{1,1}+Sl{1,1}+Sm{2,2}+Sm{3,3}     Sl{1,2}+Sm{1,2}'      Sl{1,3}+Sm{1,3}';...
     Sm{1,2}+Sl{1,2}'      2*Sm{2,2}+Sm{1,1}+Sl{2,2}+Sm{3,3}     Sl{2,3}+Sm{2,3}';...
     Sl{1,3}'+Sm{1,3}      Sl{2,3}'+Sm{2,3}     2*Sm{3,3}+Sm{1,1}+Sm{2,2}+Sl{3,3}];
  
if nargin > 7
     rho = varargin{1};
     s = varargin{2};
     Mh = massMatrix3D(rho,T,k);
%      MAT = MAT + s^2*blkdiag(Mh,Mh,Mh);  % This may be memory intensive
%      % Memory friendly approach
%      % MAT(1:end/3,1:end/3) = MAT(1:end/3,1:end/3) + s^2*Mh;
%      % MAT(end/3+1:2*(end/3),end/3+1:2*(end/3)) = MAT(end/3+1:2*(end/3),end/3+1:2*(end/3)) + s^2*Mh; 
%      % MAT(2*(end/3)+1:end,2*(end/3)+1:end) = MAT(2*(end/3)+1:end,2*(end/3)+1:end)  + s^2*Mh;  
    %Memory friendlier approach
    O=sparse(size(Mh,1),size(Mh,2));
    MAT=MAT+s^2*[Mh O O;...
                 O Mh O;...
                 O O Mh];    
end
 
[uh,~,free] = dirichletBC3D(uD,T,k);
uh=uh(:);
RHS = [fx+tx; fy+ty;fz+tz]+rhsVec-MAT*uh;

Ndof = size(fx,1);
free = [free(:); Ndof+free(:);2*Ndof+free(:)];
uh(free) = MAT(free,free)\RHS(free);
u1h = reshape(uh,[Ndof,3]);
u2h = u1h(:,2);
u3h = u1h(:,3);
u1h(:,2:3) = [];

end