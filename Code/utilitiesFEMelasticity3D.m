function [Sf,Mf,Sfd,Mfd,Dir,Free] = utilitiesFEMelasticity3D(mu,lam,T,k)
%
% [Sf,Mf,Sfd,Mfd,Dir,Free] = utilitiesFEMelasticity3D(mu,lam,T,k)
%
% Input:
% mu            : constant Lame parameter
% lam           : constant Lame parameter
% T             : Data structure. Basic FEM tetrahedrization, updated
% k             : scalar. Polynomial degree
%
% Output:
% Sf            : 3freeDOF x 3freeDOF Sparse Stiffness Matrix.
% Mf            : 3freeDOF x 3freeDOF Sparse Mass Matrix.
% Sfd           : 3freeDOF x 3dirDOF Sparse Stiffness Matrix.
% Mfd           : 3freeDOF x 3dirDOF Sparse Mass Matrix.
% Dir           : list of Dirichlet DoF in 3D system (length 3 dirDOF )
% Free          : list of Free DoF in 3D system (length 3 freeDOF )
% Last modified September 30, 2016.


rho = @(x,y,z) 0*x+1;

S = stiffnessMatricesCC3D(T,k);
S=[(2*mu+lam)*S{1,1}+mu*(S{2,2}+S{3,3})    lam*S{1,2}+mu*S{1,2}'     lam*S{1,3}+mu*S{1,3}';...
   mu*S{1,2}+lam*S{1,2}'      (2*mu+lam)*S{2,2}+mu*(S{1,1}+S{3,3})     lam*S{2,3}+mu*S{2,3}';...
   lam*S{1,3}'+mu*S{1,3}      lam*S{2,3}'+mu*S{2,3}     (2*mu+lam)*S{3,3}+mu*(S{1,1}+S{2,2})];

Mh = massMatrix3D(rho,T,k);
O=sparse(size(Mh,1),size(Mh,2));
Mass= [Mh O O;...
       O Mh O;...
       O O Mh];

dimVh=size(Mh,1);
[~,~,~,~,~,dirlist,free] = bdDOF3D(T,k);
Dir=[dirlist;dirlist+dimVh;dirlist+2*dimVh];
Free=[free;free+dimVh;free+2*dimVh];

Sf=S(Free,Free);
Mf=Mass(Free,Free);
Sfd=S(Free,Dir);
Mfd=Mass(Free,Dir);