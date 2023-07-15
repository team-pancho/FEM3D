%% C infinity exact solution for the spatial part of displacement
% Input: 
% k         : Polynomial degree. Must be >= 3.
% Output:
% rho       : Density. Vectorized function handle of x,y,z
% mu,lam    : Lame parameters. Vectorized function handle of x,y,z
% mux,lamx..: Derivatives of Lame parameters
% u,v,w     : Spatial exact solution. Vectorized function handles of x,y,z
% ux,...
% uxy,uxx,..: All first and second derivatives of u,v,w.
% Last Modified: October 23, 2018

% physical parameters
r = @(x,y,z) x.^2 + y.^2 + z.^2;
rho = @(x,y,z) 1./(1+r(x,y,z));

lam   = @(x,y,z) 1 + 1./(1+r(x,y,z));
lamx  = @(x,y,z) -2*x./(1+r(x,y,z)).^2;
lamy  = @(x,y,z) -2*y./(1+r(x,y,z)).^2;
lamz  = @(x,y,z) -2*z./(1+r(x,y,z)).^2;

mu  = @(x,y,z) 3 + cos(x.*y.*z);
mux = @(x,y,z) -y.*z.*sin(x.*y.*z);
muy = @(x,y,z) -x.*z.*sin(x.*y.*z);
muz = @(x,y,z) -x.*y.*sin(x.*y.*z);

% spacial part of exact solution
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

v     = @(x,y,z) 5*x.^2.*y.*z + 4*x.*y.^2.*z + 3*x.*y.*z.^2 + 17;
vx    = @(x,y,z) 10*x.*y.*z + 4*y.^2.*z + 3*y.*z.^2;
vy    = @(x,y,z) 5*x.^2.*z + 8*x.*y.*z + 3*x.*z.^2;
vz    = @(x,y,z) 5*x.^2.*y + 4*x.*y.^2 + 6*x.*y.*z;
vxx   = @(x,y,z) 10*y.*z;
vxy   = @(x,y,z) 10*x.*z + 8*y.*z + 3*z.^2;
vxz   = @(x,y,z) 10*x.*y + 4*y.^2 + 6*y.*z;
vyy   = @(x,y,z) 8*x.*z;
vyz   = @(x,y,z) 5*x.^2 + 8*x.*y + 6*x.*z;
vzz   = @(x,y,z) 6*x.*y;

w     = @(x,y,z) cos(2*x).*cos(3*y).*cos(z);
wx    = @(x,y,z) -2*sin(2*x).*cos(3*y).*cos(z);
wy    = @(x,y,z) -3*cos(2*x).*sin(3*y).*cos(z);
wz    = @(x,y,z) -cos(2*x).*cos(3*y).*sin(z);
wxx   = @(x,y,z) -4*w(x,y,z);
wxy   = @(x,y,z) 6*sin(2*x).*sin(3*y).*cos(z);
wxz   = @(x,y,z) 2*sin(2*x).*cos(3*y).*sin(z);
wyy   = @(x,y,z) -9*w(x,y,z);
wyz   = @(x,y,z) 3*cos(2*x).*sin(3*y).*sin(z);
wzz   = @(x,y,z) -1*w(x,y,z);
