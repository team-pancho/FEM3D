function [Uh,Vh,Wh,Sxxh,Syyh,Szzh,Syzh,Sxzh,Sxyh]= ...
    FEMViscoelasticity3D(mu,lam,rho,f,uD,sigma,T,k,time,varargin)
% [Uh,Vh,Wh,Sxxh,Syyh,Szzh,Syzh,Sxzh,Sxyh]= ...
%     FEMViscoelasticity3D(mu,lam,rho,f,uD,sigma,T,k,time)
%
% [Uh,Vh,Wh,Sxxh,Syyh,Szzh,Syzh,Sxzh,Sxyh]= ...
%     FEMViscoelasticity3D(mu,lam,rho,f,uD,sigma,T,k,time,flag)
%
% [Uh,Vh,Wh,Sxxh,Syyh,Szzh,Syzh,Sxzh,Sxyh]= ...
%     FEMViscoelasticity3D(mu,lam,rho,f,uD,sigma,T,k,time,flag,p)
%
% Input:
% mu        : First Lame parameter. Cell array {m_mu,b_mu,a_mu,nu_mu} of
%             function handles of (x,y,z)
% lam       : Second Lame parameter. Cell array {m_lam,b_lam,a_lam,nu_lam}
%             of function handles of (x,y,z)
% rho       : Mass density. Function handle of (x,y,z)
% uD        : Dirichlet BC. 3x1 cell array of function handles of (x,y,z,t)
% sigma     : Neumann BC. 3x3 cell array {Sxx,Sxy,Sxz;
%                                         Syx,Syy,Syz;
%                                         Szx,Szy,Szz}
%             of function handles of (x,y,z,t)
% f         : Force term. 3x1 cell array of function handles of (x,y,z,t)
% T         : Enhanced tetrahedrization
% k         : Polynomial degree
% time      : Equipartitioned temporal grid of [0,T] (vector of 1 x Nt)
% flag      : Optional. 0 (Default) for stress at nodes, 1 for ave. stress
% p         : Optional. Time stepping method. Function handle for the CQ 
%             method or matrix for RKCQ. Default is trapezoidal rule
%
% Output:
% Uh,Vh,Wh  : Ndof x Nt array of FE cofficients for displacement in TD.
% Sxxh,...  : Components of stress in TD. dimP_{k-1}*Nelt x Nt or Nelt x Nt
%             array depending on whether flag is 0 or 1
%
% Last modified: Oct 24, 2018

nofarg = 9; % number of non optional arguments
% If 1st optional argument exists and = 1, then average stress is computed.
avestressFlag = 0;
if nargin >= nofarg+1
    avestressFlag = varargin{1};
end
% Handling 2nd optional argument for CQ/RKCQ
p = @(z) 2*(1-z)./(1+z); % Trapezoidal rule is default method
if nargin>=nofarg+2
    p = varargin{2};
end
% Determining the time stepping method
if isa(p,'function_handle')
    multistepMethod = true;
elseif isa(p,'numeric')
    A = p; % Runge-Kutta matrix
    multistepMethod = false;
    S = size(A,1); % number of stages
    c = sum(A,2);
end
% Dimensions
Nelt = size(T.elements,2);
Ndof = dimFEMspace(T,k);
km1 = k-1;
dimPkm1Nelt = (km1+3)*(km1+2)*(km1+1)/6*size(T.elements,2);

% Time-step
dt = time(2) - time(1);

% Sampling the data for CQ or RKCQ methods
if multistepMethod % implementing CQ
    Nt = length(time);
    % Generate input data for transfer function
    Y = zeros(3*Ndof,Nt); Z = zeros(3*Ndof,Nt);
    parfor n=1:Nt
        t=time(n);
        [fh,trh,uDh]=sampleData(f,sigma,uD,T,k,t)
        Y(:,n) = uDh;
        Z(:,n) = fh + trh;
    end
    X = [Y;Z];
else % implementing RKCQ
    Nt = length(time) - 1;
    Y = zeros(3*Ndof,S,Nt); Z = zeros(3*Ndof,S,Nt);
    X = zeros(6*Ndof,S,Nt);
    parfor n=1:Nt
        tn=time(n);
        for ns=1:S
            t = tn + dt * c(ns); %#ok
            [fh,trh,uDh]=sampleData(f,sigma,uD,T,k,t);
            Z(:,ns,n) = fh+trh;
            Y(:,ns,n) = uDh;
        end
    end
    % Reshape the data
    for ns=1:S
        X(:,ns,:) = [Y(:,ns,:);Z(:,ns,:)];
    end
    X = reshape(X,[S*6*Ndof,Nt]);
end

% Lame parameters are unpacked
m_mu = mu{1}; b_mu = mu{2}; a_mu = mu{3}; nu_mu = mu{4};
m_lam = lam{1}; b_lam = lam{2}; a_lam = lam{3}; nu_lam = lam{4};
muts = @(x,y,z,s) (m_mu(x,y,z) + ...
    s.^nu_mu(x,y,z).*b_mu(x,y,z))./(1 + s.^nu_mu(x,y,z).*a_mu(x,y,z));
lamts = @(x,y,z,s) (m_lam(x,y,z) + ...
    s.^nu_lam(x,y,z).*b_lam(x,y,z))./(1 + s.^nu_lam(x,y,z).*a_lam(x,y,z));
% Construct transfer function for time-stepping
Transfer = @(s,X) TFFEMelasticity3D(muts,lamts,rho,X,T,k,s,avestressFlag);

if avestressFlag
    dimStress = Nelt;
else
    dimStress = dimPkm1Nelt;
end
% CQ outputs all the information in one long vector
if multistepMethod
    allOut=CQoperator(Transfer,X,dt,3*Ndof+6*dimStress,p);
else % Runge-Kutta method is being implemented
    allOut=RKCQoperator(Transfer,X,dt,3*Ndof+6*dimStress,0,p);
end
% Indices to extract displacement and stress
block1=@(x) (1+(x-1)*Ndof):(x*Ndof);
block2=@(x) (1+(x-1)*dimPkm1Nelt):(x*dimStress);
% Unpacking displacement and stess
Uh = allOut(block1(1),:);
Vh = allOut(block1(2),:);
Wh = allOut(block1(3),:);
Sxxh = allOut(Ndof*3+block2(1),:);
Syyh = allOut(Ndof*3+block2(2),:);
Szzh = allOut(Ndof*3+block2(3),:);
Syzh = allOut(Ndof*3+block2(4),:);
Sxzh = allOut(Ndof*3+block2(5),:);
Sxyh = allOut(Ndof*3+block2(6),:);
end

function [fh,trh,uDh]=sampleData(f,sigma,uD,T,k,t)
% This is a helper function. Given time value t, it returns the sampled
% data for force, traction and dirichlet tests

% Load function
ft = cell(1,3);
for i=1:3
    ft{i} = @(x,y,z) f{i}(x,y,z,t);
end
fh = loadVector3D(ft,T,k);
fh = fh(:);

% Tranction
sigmat = cell(3,3);
for i=1:3
    for j=1:3
        sigmat{i,j} = @(x,y,z) sigma{i,j} (x,y,z,t);
    end
end
trh=neumannBC3D(sigmat,T,k);
trh = trh(:);

% Dirichlet condition
uDt = cell(1,3);
for i=1:3
    uDt{i} = @(x,y,z) uD{i} (x,y,z,t);
end
[uDh,~,~] = dirichletBC3D(uDt,T,k);
uDh = uDh(:);
end