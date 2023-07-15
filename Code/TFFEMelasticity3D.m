function sol = TFFEMelasticity3D(muts,lamts,rho,Data,T,k,s,varargin)
% [Uhx;Uhy;Uhz;Sigxx;Sigyy;Sigzz;Sigyz;Sigxz;Sigxy] ...
%           = TFFEMelasticity3D(muts,lamts,rho,[uDh;rhs],T,k,s,flag)
%                                          
% Input:
% muts, lamts   : Viscoelastic Lame parameters. Functions of (x,y,z,s)
% rho           : Mass density. Function of (x,y,z)
% [uDh;rhs]     : 6*Ndof x 1. Dirichlet BC and RHS (force term + Neumann BC)
% T             : Enhanced tetrahedrization
% k             : Polynomial degree
% s             : Complex number. Frequency parameter in Laplace domain
% flag          : 0 for stress at nodes, 1 for avereged stress at each
%                 element
%
% Output:
% [Uhx;...;Sigxx;...] : (3*Ndof + 6*dimStress) x 1 vector where
%                        Uhx,... : Ndof x 1 vectors for FE coefficients of 
%                           the displacement. 
%                        Sigxx,...: Coefficients of stress in piecewise 
%                           P_{k-1} space or averaged stress on each element
%                           where dimStress is either dimP_{k-1}*Nelt or 
%                           Nelt depending on whether flag is 0 or 1.
%
% Last modified : Oct 24, 2018

% Extracting data
Ndof = dimFEMspace(T,k);
uh = Data(1:3*Ndof);
rhs = Data(3*Ndof+1:6*Ndof);
% Constructing stiffness and mass matrices
mu = @(x,y,z) muts(x,y,z,s);
lam = @(x,y,z) lamts(x,y,z,s);
Sm = stiffnessMatrices3D(mu,T,k);
Sl = stiffnessMatrices3D(lam,T,k);
Mh = massMatrix3D(rho,T,k);

Sh = [2*Sm{1,1}+Sl{1,1}+Sm{2,2}+Sm{3,3}     Sl{1,2}+Sm{1,2}.'      Sl{1,3}+Sm{1,3}.';...
      Sm{1,2}+Sl{1,2}.'      2*Sm{2,2}+Sm{1,1}+Sl{2,2}+Sm{3,3}     Sl{2,3}+Sm{2,3}.';...
      Sl{1,3}.'+Sm{1,3}      Sl{2,3}.'+Sm{2,3}     2*Sm{3,3}+Sm{1,1}+Sm{2,2}+Sl{3,3}];

O = sparse(size(Mh,1),size(Mh,2));
Mh = [Mh O O;...
      O Mh O;...
      O O Mh];
% Getting Dricihlet and free indices
[~,~,~,~,~,dir,free]= bdDOF3D(T,k);
dir  = [dir(:) ; Ndof+dir(:) ; 2*Ndof+dir(:)];
free = [free(:); Ndof+free(:); 2*Ndof+free(:)]; 
% Solving for the displacement 
MSh = s^2*Mh + Sh;
uh(free) = MSh(free,free)\...
    (rhs(free) - MSh(free,dir)*uh(dir));

u1h = reshape(uh,[Ndof,3]);
u2h = u1h(:,2);
u3h = u1h(:,3);
u1h(:,2:3) = [];
% In case of the optional argument is passed in and is equal to 1 then 
% average stress is computed
nofarg = 7; % Number of non-optional arguments
if nargin >= nofarg+1
    avestress_flag = varargin{1};
else
    avestress_flag = 0;
end
if avestress_flag % compute average stress
    [sigxxh,sigyyh,sigzzh,sigyzh,sigxzh,sigxyh] = ...
        stressPostprocessing(u1h,u2h,u3h,lam,mu,T,k,1);
    % Combining displacement and stress as a single vector output
    sol = [u1h;u2h;u3h;sigxxh.';sigyyh.';sigzzh.';sigyzh.';sigxzh.';sigxyh.'];
    return
else % compute stress at nodes
    [sigxxh,sigyyh,sigzzh,sigyzh,sigxzh,sigxyh] = ...
        stressPostprocessing(u1h,u2h,u3h,lam,mu,T,k);
end
% Combining displacement and stress as a single vector output
sol = [u1h;u2h;u3h;sigxxh;sigyyh;sigzzh;sigyzh;sigxzh;sigxyh];
end
