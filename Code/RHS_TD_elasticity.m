function [RHSfree] = RHS_TD_elasticity(t,U,T,k,Sf,Sfd,Mfd,g0)
%
% RHS_TD_elasticity(t,U,T,k,Sf,Sfd,Mfd,g0)
%
% Input:
% t             : time (one instance)
% U             : vector of coefficients of uh from previous timestep
% T             : Data structure. Basic FEM tetrahedrization, updated
% k             : scalar. Polynomial degree
% Sf            : 3freeDOF x 3freeDOF sparse stiffness matrix
% Sfd           : 3freeDOF x 3dirDOF Sparse Stiffness Matrix.
% Mfd           : 3freeDOF x 3dirDOF Sparse Mass Matrix.
% g0            : array of six functions {uxD,uyD,uzD,uxD_tt,uyD_tt,uzD_tt}
%
% Output:
% Yfree         : Values for free DOF on RHS of TD elasticity.
% Last modified July 7, 2015.

d3=size(U,1);
for i=1:6
   f{i}= @(x,y,z) g0{i}(t,x,y,z);
end
[Udir,dir]=dirichletBC3D({f{1} f{2} f{3}},T,k);
Udir=Udir(dir,:);
Udir=Udir(:);

Udirtt=dirichletBC3D({f{4} f{5} f{6}},T,k);
Udirtt=Udirtt(dir,:);
Udirtt=Udirtt(:);

RHSfree=Sf*U-Mfd*Udir-Sfd*Udirtt;


