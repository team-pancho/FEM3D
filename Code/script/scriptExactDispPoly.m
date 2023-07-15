%% Polynomial exact solution for the spatial part of displacement
% Input: 
% k         : Polynomial degree. Must be >= 3.
% Output:
% rho       : Density. Degree 1 vectorized function handle of x,y,z.
% mu,lam    : Lame parameters. Degree 1 vect. function handle of x,y,z.
% mux,lamx..: Derivatives of Lame parameters
% u,v,w     : Spatial exact solution. Degree k-1 function handles of x,y,z.
% ux,...
% uxy,uxx,..: All first and second derivatives of u,v,w.
% Note that this setup will give degree k-1 stress, degree k-2 forcing term
% Last Modified: October 23, 2018

% physical parameters
assert(k>=3,'Polynomial degree should be >=3');
% Well posed when mesh has x,y,z>=0
% symbolic calculation
syms x y z

rho = x + 2*y + 10*z + 7;
rho = matlabFunction(rho,'Vars',[x y z]);
lam = x + sqrt(11)*y + 4*z + 9;

lamx = matlabFunction(diff(lam,x),'Vars',[x y z]);
lamy = matlabFunction(diff(lam,y),'Vars',[x y z]);
lamz = matlabFunction(diff(lam,z),'Vars',[x y z]);
lam = matlabFunction(lam,'Vars',[x y z]);

mu  = 2*y + 5*x + 7*z + 5;

mux = matlabFunction(diff(mu,x),'Vars',[x y z]);
muy = matlabFunction(diff(mu,y),'Vars',[x y z]);
muz = matlabFunction(diff(mu,z),'Vars',[x y z]);
mu = matlabFunction(mu,'Vars',[x y z]);

u = 5*x.^(k-1) - 2*y^(k-1) + x.*z;
v = sqrt(5)*y.^(k-1) + 10*z^(k-1) - x.*y;
w = -4*x.^(k-1) + 6*z^(k-1) + x + y*z - z^2;

ux = matlabFunction(diff(u,x),'Vars',[x y z]);
uy = matlabFunction(diff(u,y),'Vars',[x y z]);
uz = matlabFunction(diff(u,z),'Vars',[x y z]);
uxx = matlabFunction(diff(ux,x),'Vars',[x y z]);
uxy = matlabFunction(diff(ux,y),'Vars',[x y z]);
uxz = matlabFunction(diff(ux,z),'Vars',[x y z]);
uyy = matlabFunction(diff(uy,y),'Vars',[x y z]);
uyz = matlabFunction(diff(uy,z),'Vars',[x y z]);
uzz = matlabFunction(diff(uz,z),'Vars',[x y z]);

vx = matlabFunction(diff(v,x),'Vars',[x y z]);
vy = matlabFunction(diff(v,y),'Vars',[x y z]);
vz = matlabFunction(diff(v,z),'Vars',[x y z]);
vxx = matlabFunction(diff(vx,x),'Vars',[x y z]);
vxy = matlabFunction(diff(vx,y),'Vars',[x y z]);
vxz = matlabFunction(diff(vx,z),'Vars',[x y z]);
vyy = matlabFunction(diff(vy,y),'Vars',[x y z]);
vyz = matlabFunction(diff(vy,z),'Vars',[x y z]);
vzz = matlabFunction(diff(vz,z),'Vars',[x y z]);

wx = matlabFunction(diff(w,x),'Vars',[x y z]);
wy = matlabFunction(diff(w,y),'Vars',[x y z]);
wz = matlabFunction(diff(w,z),'Vars',[x y z]);
wxx = matlabFunction(diff(wx,x),'Vars',[x y z]);
wxy = matlabFunction(diff(wx,y),'Vars',[x y z]);
wxz = matlabFunction(diff(wx,z),'Vars',[x y z]);
wyy = matlabFunction(diff(wy,y),'Vars',[x y z]);
wyz = matlabFunction(diff(wy,z),'Vars',[x y z]);
wzz = matlabFunction(diff(wz,z),'Vars',[x y z]);

u = matlabFunction(u,'Vars',[x y z]);
v = matlabFunction(v,'Vars',[x y z]);
w = matlabFunction(w,'Vars',[x y z]);