%% scriptExactSolutionFracVisco.m
% A script to generate an exact solution for fractional isotropic
% generalized Zener model for viscoelastic wave propagation problem.
% Different viscoelastic models can be tested via different input.
%
% Input:
% a,b,m         : Viscoelastic parameters. Scalars such that b > a*m
% nu            : Fractional power. Scalar in [0,1]
% deg           : Degree of H(t). Integer >= 3
% Spatial func.s: Should be defined before calling the script
%                 Mass density : rho
%                 Lame and derv: mu, lam, mux, muy, muz, lamx, lamy, lamz
%                 Displacement : u, v, w,
%                 Derv. of disp: ux, uy, uz, uxx, uxy, uxx, uyy, uyz, uzz,
%                                and same for v, w
%
% Output:
% U,V,W     : Seperable exact solution. Function handles of x,y,z,t
%             U = u(x)h(t), V=v(x)h(t), W=w(x)h(t)
% Ux,Uy,... : All derivatives of U,V,W
% Sxx,...   : Symmetric stress tensor. 6 function handles of x,y,z,t.
% Fx,Fy,Fz  : Forcing term. Function handles of x,y,z,t
% Mu, Lam   : Lame parameters. 4 x 1 cell array of function handles of
%             (x,y,z,t). See FEMViscoelasticity3D's input
%
% Last modified: February 1, 2020

%% fractional Zener model and temporal solution parameters
% These parameters help to obtain different models
% c1 - a*c0 > 0, c1=b, c0=m, hence b - a*m > 0
% Example input:  a = 1; b = 2.5; m = 1; nu = 0.5; deg = 3;
if ~exist('external','var')
    a = input('Viscoparameter a: ');
    b = input('Viscoparameter b: ');
    m = input('Viscoparameter m (so that b>a*m): ');
    nu = input('Fractional power nu: ');
    deg = input('Degree of temporal solution: ');
end
%% Derived quantities: Spatial part of forcing term and stress
div=@(x,y,z) ux(x,y,z)+vy(x,y,z)+wz(x,y,z);
divx=@(x,y,z) uxx(x,y,z)+vxy(x,y,z)+ wxz(x,y,z);
divy=@(x,y,z) uxy(x,y,z)+vyy(x,y,z)+wyz(x,y,z);
divz=@(x,y,z) uxz(x,y,z)+vyz(x,y,z)+wzz(x,y,z);

sxx=@(x,y,z) 2*mu(x,y,z).*ux(x,y,z)+lam(x,y,z).*div(x,y,z);
syy=@(x,y,z) 2*mu(x,y,z).*vy(x,y,z)+lam(x,y,z).*div(x,y,z);
szz=@(x,y,z) 2*mu(x,y,z).*wz(x,y,z)+lam(x,y,z).*div(x,y,z);
sxy=@(x,y,z) mu(x,y,z).*(vx(x,y,z)+uy(x,y,z));
sxz=@(x,y,z) mu(x,y,z).*(wx(x,y,z)+uz(x,y,z));
syz=@(x,y,z) mu(x,y,z).*(wy(x,y,z)+vz(x,y,z));
fx=@(x,y,z) -(lam(x,y,z).*divx(x,y,z)+lamx(x,y,z).*div(x,y,z)...
    +2*mu(x,y,z).*uxx(x,y,z)+2*mux(x,y,z).*ux(x,y,z)...
    +mu(x,y,z).*(vxy(x,y,z)+uyy(x,y,z))+muy(x,y,z).*(vx(x,y,z)+uy(x,y,z))...
    +mu(x,y,z).*(wxz(x,y,z)+uzz(x,y,z))+muz(x,y,z).*(wx(x,y,z)+uz(x,y,z)));
fy=@(x,y,z) -(lam(x,y,z).*divy(x,y,z)+lamy(x,y,z).*div(x,y,z)...
    +2*mu(x,y,z).*vyy(x,y,z)+2*muy(x,y,z).*vy(x,y,z)...
    +mu(x,y,z).*(wyz(x,y,z)+vzz(x,y,z))+muz(x,y,z).*(wy(x,y,z)+vz(x,y,z))...
    +mu(x,y,z).*(uxy(x,y,z)+vxx(x,y,z))+mux(x,y,z).*(uy(x,y,z)+vx(x,y,z)));
fz=@(x,y,z) -(lam(x,y,z).*divz(x,y,z)+lamz(x,y,z).*div(x,y,z)...
    +2*mu(x,y,z).*wzz(x,y,z)+2*muz(x,y,z).*wz(x,y,z)...
    +mu(x,y,z).*(uxz(x,y,z)+wxx(x,y,z))+mux(x,y,z).*(uz(x,y,z)+wx(x,y,z))...
    +mu(x,y,z).*(vyz(x,y,z)+wyy(x,y,z))+muy(x,y,z).*(vz(x,y,z)+wy(x,y,z)));

%% Temporal part of exact solution
h = @(t) t.^deg;
hp = @(t) deg*t.^(deg-1);
hpp = @(t) deg*(deg-1)*t.^(deg-2);
prec_mlf = 10; % accuracy of the computation of mlf. 10 should be enough.
if a == 0
    sigmaT = @(t) m*t.^(deg) + b*gamma(deg+1)/gamma(deg+1-nu)*t.^(deg-nu);
else
    sigmaT = @(t) m*gamma(deg+1)/a*t.^(nu+deg).*mlf(nu,nu+deg+1,-t.^nu/a,prec_mlf)...
        +b*gamma(deg+1)/a*t.^deg.*mlf(nu,deg+1,-t.^nu/a,prec_mlf);
end

%% Combining spatial and temporal solutions
% exact stress
Sxx = @(x,y,z,t) sxx(x,y,z).*sigmaT(t);
Syy = @(x,y,z,t) syy(x,y,z).*sigmaT(t);
Szz = @(x,y,z,t) szz(x,y,z).*sigmaT(t);
Sxy = @(x,y,z,t) sxy(x,y,z).*sigmaT(t);
Sxz = @(x,y,z,t) sxz(x,y,z).*sigmaT(t);
Syz = @(x,y,z,t) syz(x,y,z).*sigmaT(t);

% exact forcing term
Fx = @(x,y,z,t) rho(x,y,z)*hpp(t).*u(x,y,z) + fx(x,y,z).*sigmaT(t);
Fy = @(x,y,z,t) rho(x,y,z)*hpp(t).*v(x,y,z) + fy(x,y,z).*sigmaT(t);
Fz = @(x,y,z,t) rho(x,y,z)*hpp(t).*w(x,y,z) + fz(x,y,z).*sigmaT(t);
% exact displacement
U = @(x,y,z,t) u(x,y,z)*h(t);
V = @(x,y,z,t) v(x,y,z)*h(t);
W = @(x,y,z,t) w(x,y,z)*h(t);

Ux = @(x,y,z,t) ux(x,y,z)*h(t);
Vx = @(x,y,z,t) vx(x,y,z)*h(t);
Wx = @(x,y,z,t) wx(x,y,z)*h(t);

Uy = @(x,y,z,t) uy(x,y,z)*h(t);
Vy = @(x,y,z,t) vy(x,y,z)*h(t);
Wy = @(x,y,z,t) wy(x,y,z)*h(t);

Uz = @(x,y,z,t) uz(x,y,z)*h(t);
Vz = @(x,y,z,t) vz(x,y,z)*h(t);
Wz = @(x,y,z,t) wz(x,y,z)*h(t);

Mu = cell(1,4);
Lam = cell(1,4);
Mu{1} = @ (x,y,z) m*mu(x,y,z);
Mu{2} = @ (x,y,z) b*mu(x,y,z);
Mu{3} = @ (x,y,z) a + 0*x;
Mu{4} = @ (x,y,z) nu + 0*x;
Lam{1} = @ (x,y,z) m*lam(x,y,z);
Lam{2} = @ (x,y,z) b*lam(x,y,z);
Lam{3} = @ (x,y,z) a + 0*x;
Lam{4} = @ (x,y,z) nu + 0*x;