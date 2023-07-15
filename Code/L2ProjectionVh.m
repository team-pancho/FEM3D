function projVh = L2ProjectionVh(f,T,k)

% projVh = L2ProjectionVh(f,T,k);
% 
% Input :
%    f  : vectorized function of three variables
%    T  : enhanced triangulation
%    k  : polynomial degree
%
% Output:
%    projVh : L2 projection of f(x,y,z) onto the finite element space
%
% Last Modified: May 1, 2015

c=@(x,y,z) 1+0*x;
Mh = massMatrix3D(c,T,k);
RHS = loadVector3D(f,T,k);

projVh = Mh\RHS;

return