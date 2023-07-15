function [u1h,u2h,u3h,sxxh,syyh,szzh,syzh,sxzh,sxyh]...
                 = FEMelasticity3D(mu,lam,rho,f,uD,sigma,T,k,s)
% [u1h,u2h,u3h,sxxh,syyh,szzh,syzh,sxzh,sxyh]...
%                  = FEMelasticity3D(mu,lam,rho,{fx,fy,fz},{u,v,w},...
%                                    {sxx,sxy,sxz,syy,syz,szz},T,k,s)
% Input:
% mu,lam, rho   : Vectorized function handles of three variables
% fx, fy, fz    : Vectorized functions of three variables (load vector)
% ux, uy, uz    : Vectorized functions of three variables (Dirichlet BC)
% sigxx,...     : Six vectorized functions of three variables
% T             : Data structure. Basic FEM triangulation
% k             : scalar. Polynomial degree
% s             : Complex wave number (s=-1i k)
%
% Output:
% u1h,u2h,u3h   : Ndof x 1 Vectors. Vh coefficients of the displacement
%                 components.
% sxxh,...      : dim P_{k-1} times N_elts vectors.  Coefficients of the DG
%                 stress post-processed using the FEM displacement solution
%
% Last modified: November 21, 2017

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

Sm = stiffnessMatrices3D(mu,T,k);
Sl = stiffnessMatrices3D(lam,T,k);

MAT=[2*Sm{1,1}+Sl{1,1}+Sm{2,2}+Sm{3,3}  Sl{1,2}+Sm{1,2}.'   Sl{1,3}+Sm{1,3}.';...
     Sm{1,2}+Sl{1,2}.'  2*Sm{2,2}+Sm{1,1}+Sl{2,2}+Sm{3,3}   Sl{2,3}+Sm{2,3}.';...
     Sl{1,3}.'+Sm{1,3}  Sl{2,3}.'+Sm{2,3}  2*Sm{3,3}+Sm{1,1}+Sm{2,2}+Sl{3,3}];
  
Mh = massMatrix3D(rho,T,k);
O=sparse(size(Mh,1),size(Mh,2));
MAT=MAT+s^2*[Mh O O;...
             O Mh O;...
             O O Mh];    
 
[uh,~,free] = dirichletBC3D(uD,T,k);
uh=uh(:);
RHS = [fx+tx;fy+ty;fz+tz]-MAT*uh;

Ndof = dimFEMspace(T,k);
free = [free(:); Ndof+free(:);2*Ndof+free(:)];
uh(free) = MAT(free,free)\RHS(free);
u1h = reshape(uh,[Ndof,3]);
u2h = u1h(:,2);
u3h = u1h(:,3);
u1h(:,2:3) = [];

[sxxh,syyh,szzh,syzh,sxzh,sxyh]=stressPostprocessing(u1h,u2h,u3h,lam,mu,T,k);

end