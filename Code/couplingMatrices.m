function [MAT,b] = couplingMatrices(T,Tp,Tm,ep,em,rho,lam,mu,k,s,g,G,vinc)

% [Mh,b]= couplingMatrices(T,Tp,Tm,ep,em,rho,lam,mu,k,s,g,G,vinc)
% Input:
%   Tp, Tm          : Data structure. Enhanced FEM triangulation
%   ep, em          : Embeding operators
%   mu,lam          : Vect function of 3 vars (elastic parameters)
%   rho             : Vect function of 3 vars (density)
%   g               : 1x3 cellarray w vect fns of 3 vars (jump in velocity)
%   G               : 1x6 cellarray w vect fnd of 3 vars (jump in stress)
%                         Order: xx,xy,xz,yy,yz,zz
%   k               : Polynomial degree
%   s               : Complex frequency in Laplace domain
%   vinc            : 1x4 cell array w vect fns of 3 vers (inc wave & grad)
%
% Output:
%   Mh,b            :
%
% Last modified : November 18, 2016

% Elasticity in Tm

Sm = stiffnessMatrices3D(mu,Tm,k);
Sl = stiffnessMatrices3D(lam,Tm,k);
Mel = massMatrix3D(rho,Tm,k);
Sel = [2*Sm{1,1}+Sl{1,1}+Sm{2,2}+Sm{3,3}, Sl{1,2}+Sm{1,2}.',...
                                                    Sl{1,3}+Sm{1,3}.';...
        Sm{1,2}+Sl{1,2}.', 2*Sm{2,2}+Sm{1,1}+Sl{2,2}+Sm{3,3},...
                                                    Sl{2,3}+Sm{2,3}.';...
        Sl{1,3}.'+Sm{1,3}, Sl{2,3}.'+Sm{2,3},...
                                    2*Sm{3,3}+Sm{1,1}+Sm{2,2}+Sl{3,3}];
O = sparse(size(Mel,1),size(Mel,2));
Mel = [Mel O O;...
       O Mel O;...
       O O Mel];
EL = s^2*Mel + Sel;

% Acoustics in Tp

unit = @(x,y,z) 0*x+1;
Mac = massMatrix3D(unit,Tp,k);
Sac = stiffnessMatricesCC3D(Tp,k);
AC = s^2*Mac + Sac{1,1}+Sac{2,2}+Sac{3,3};
[~,MBd] = neumannBC3D(unit,Tp,k,unit);  % b.form impedance condition
AC = AC+s*MBd;

% Coupling terms 

d = dimFEMspace(T,k);
dm = dimFEMspace(Tm,k);
dp = dimFEMspace(Tp,k);
I = interfaceMatrices3D(Tm,ep,em,[d,dp,dm],k,g,G);
C = [I{2};I{3};I{4}];

MAT = [EL   s*C;
       -s*C'  AC];

% Right-hand side

bext = s*neumannBC3D(vinc{1},Tp,k)+neumannBC3D({vinc{2:4}},Tp,k);
b = [I{5}; -I{6} + bext];


end