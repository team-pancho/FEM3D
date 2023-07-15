% BENCHMARK SCRIPT
% linear elasticity, nonhomogeneous isotropic,
% exact solution known
% Last Modified: November 21, 2016

clear;
close all;

scriptExactSolutionElasticity
scriptDerivedQuantitiesElasticity

k=input('Polynomial degree: ');
DBC = input('Dirichlet faces (vector):');
for j=1:5
    T = meshCube(j,j,j,DBC);
    [uh,vh,wh,sxxh,syyh,szzh,syzh,sxzh,sxyh]=FEMelasticity3D(mu,lam,rho,...
        {fx,fy,fz},{u,v,w},{sxx,sxy,sxz,syy,syz,szz},T,k,s);
    [ehu(j),hhu(j)]=errorFEM3D({u,ux,uy,uz},uh,T,k);
    [ehv(j),hhv(j)]=errorFEM3D({v,vx,vy,vz},vh,T,k);
    [ehw(j),hhw(j)]=errorFEM3D({w,wx,wy,wz},wh,T,k);
    esxxh(j) = errorDFEM3D(sxx,sxxh,T,k);
    esyyh(j) = errorDFEM3D(syy,syyh,T,k);
    eszzh(j) = errorDFEM3D(szz,szzh,T,k);
    esyzh(j) = errorDFEM3D(syz,syzh,T,k);
    esxzh(j) = errorDFEM3D(sxz,sxzh,T,k);
    esxyh(j) = errorDFEM3D(sxy,sxyh,T,k);    
    h(j) = 1/j;
end

eh = ehu+ehv+ehw;
hh = hhu+hhv+hhw;
esigma = esxxh+esyyh+eszzh+esyzh+esxzh+esxyh;
% H1 errors
loglog(h,hh, '-o');
hold on
loglog(h,hh(1)*h.^k,'-');

% L2 errors
loglog(h, eh, '-or',h, esigma,'-p');
hold on
loglog(h,eh(1)*h.^(k+1),'-r');

legend('u_h H^1 error', sprintf('%s %s','Order',num2str(k)),'u_h L^2 error',...
       '\sigma_h L^2 error',sprintf('%s %s','Order',num2str(k+1)),...
       'Location', 'SouthEast')
 