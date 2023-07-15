function I = interfaceMatrices3D(Tm, ep, em, dim, k, varargin)

% I = interfaceMatrix3D(Tm, ep, em, [d,dp,dm], k)
% I = interfaceMatrix3D(Tm, ep, em, [d,dp,dm], k, g , G)
% 
% Inputs:
%       Tm : enhanced mesh of interior domain
%       ep : vector embedding Tp into T
%       em : vector embedding Tm into T
%       [d,dp,dm]: dimensions of Pk FEM spaces Vh, Vhp, Vhm
%       k  : polynomial degree
%       g  : 1x3 cell array with 3 vect fns of 3 variables
%       G  : 1x6 cell array with 6 vect fns of 3 variables
% Output:
%       I  : Cell Array of four interface matrices
%               \int_Gint P_i^- P_j^+ \tau \, \tau\in {1,nx,ny,nz}
%            ... or two added interface vectors
%               \int_Gint (G n) P_i^- 
%               \int_Gint (g.n) P_i^+
% 
% Last Modified: October 28, 2016

d = dim(1);
dp = dim(2);
dm = dim(3);

% copy on Tm.interface onto Tm.neumann
Tm = interface(Tm);
Tm.faces(4,Tm.faces(4,:)==2) = 5;
Tm.faces(4,Tm.faces(4,:)==3) = 2;
Tm.oldneumann = Tm.neumann;
Tm.neumann = Tm.interface;

onefunc = @(x,y,z) 1 + 0*x;
zerofunc = @(x,y,z) 0*x;

I{1} = sparse(dm,dp); I{2} = I{1}; I{3} = I{1}; I{4} = I{1};
[~,A] = neumannBC3D(zerofunc,Tm,k,onefunc);
B = sparse(dm,d);   
B(:,em) = A;
I{1} = B(:,ep);

[~,A] = neumannBC3D(zerofunc,Tm,k,{onefunc,zerofunc,zerofunc});
B = sparse(dm,d);   
B(:,em) = A;
I{2} = B(:,ep);

[~,A] = neumannBC3D(zerofunc,Tm,k,{zerofunc,onefunc,zerofunc});
B = sparse(dm,d);   
B(:,em) = A;
I{3} = B(:,ep);

[~,A] = neumannBC3D(zerofunc,Tm,k,{zerofunc,zerofunc,onefunc});
B = sparse(dm,d);   
B(:,em) = A;
I{4} = B(:,ep);

if nargin==7
    bG = neumannBC3D(varargin{2},Tm,k); I{5}=bG(:);
    bpg=zeros(d,1);
    bpg(em)=neumannBC3D(varargin{1},Tm,k);
    I{6} = bpg(ep);
end

