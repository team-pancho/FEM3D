% linear piezo-elasticity, nonhomogeneous isotropic,
% exact solution known
% we assume the third order piezoelectric and 
% second order dielectric tensors are constant for now
% Last Modified: December 7, 2017

clear; clc; close all;

% input parameters and exact solution
r = @(x,y,z) x.^2 + y.^2 + z.^2;
zz  = @(x,y,z) 0*x;
steady_trig %exact solutions are in the Control folder
derivedSteadyPiezo %also in the Control folder

% discrete parameters
k=input('Polynomial degree: ');
DBC=input('Dirichlet: ');
s = input('s (1 or 0): ');
% switch BCs
%     case 1
%         DBC=1:6;
%     case 2
%         DBC=[];
%     case 3
%         DBC=[1 2 4];
% end


if s
  fx = @(x,y,z) fx(x,y,z) + rho(x,y,z).*u(x,y,z);
  fy = @(x,y,z) fy(x,y,z) + rho(x,y,z).*v(x,y,z);
  fz = @(x,y,z) fz(x,y,z) + rho(x,y,z).*w(x,y,z);
end


for j=1:5
    T = meshCube(j,j,j,[DBC]);
    [uh,vh,wh,psih]=FEMpiezoelasticity3D(mu,lam,rho,e,kap,...
        {fx,fy,fz,fpsi},{u,v,w,psi},{sxx,sxy,sxz,syy,syz,szz,etax,etay,...
         etaz},T,k,s);
    [ehu(j),hhu(j)]=errorFEM3D({u,ux,uy,uz},uh,T,k);
    [ehv(j),hhv(j)]=errorFEM3D({v,vx,vy,vz},vh,T,k);
    [ehw(j),hhw(j)]=errorFEM3D({w,wx,wy,wz},wh,T,k);
    [ehpsi(j), hhpsi(j)] = errorFEM3D({psi,psix,psiy,psiz},psih,T,k);
    h(j) = 1/j;
end

eh = ehu+ehv+ehw;
hh = hhu+hhv+hhw;

% H1 errors
loglog(h,hh, '-o');
hold on
loglog(h,hhpsi,'-x');
loglog(h,hh(1)*h.^k,'-');


% L2 errors
loglog(h, eh, '-or');
hold on
loglog(h, ehpsi,'-xr');
loglog(h,eh(1)*h.^(k+1),'-r');

legend('u H1 error','psi H1 error', 'Order k',...
       'u L2 error','psi L2 error','Order k+1','Location','SouthEast')