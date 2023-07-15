function [U0,V0]=ICelasticity(u0,v0,T,k)

% [U0,V0]=ICelasticity({u0x,u0y,u0z},{v0x,v0y,v0z},T,k)
%
% Input: 
%     {u0x,u0y,u0z} : vectorized functions (initial conditions)
%     {v0x,v0y,v0z} : vectorized functions (initial conditions)
%     T             : basic or expanded tetrahedrization data structure
%     k             : polynomial degree
% Output:
%     U0            : coefficients in k-dim FEM space of {u0} by L2 proj
%     V0            : coefficients in k-dim FEM space of {v0} by L2 proj
% Last Modified: July 16th, 2015

c=@(x,y,z) 1+0*x;
Mh = massMatrix3D(c,T,k);

u0v0rhs=loadvector3D({u0{1},u0{2},u0{3},v0{1},v0{2},v0{3}},T,k);

projVh = Mh\u0v0rhs;

U0=projVh(:,1:3);  U0=U0(:);
V0=projVh(:,4:6);  V0=V0(:);
return