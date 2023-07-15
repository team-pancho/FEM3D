function [gh,MBd]=neumannBC3D(g,T,k,varargin)

% [gh]=neumannBC3D(g,T,k)
% [gh,MBd]=neumannBC3D(g,T,k,lambda)
% Input:
%     g      : a single vectorized function of three variables
%                       or
%              array of three vectorized functions of three vars {gx,gy,gz}
%                       or
%              array of six vectorized function of three variables
%                       {sigmaxx,sigmaxy,sigmaxz,sigmayy,sigmayz,sigmazz}
%              3x3 cell array of vectorized functions of three variables 
%     T      : enhanced triangulation
%     k      : polynomial degree
%     lambda : array of vectorized functions of three variables 
%                       or
%              a single vectorized function of three variables
%               (optional)
% Output:
%     gh     : dimVh x 1boundary vector int_{Gamma}(g or g dot n) \phi_i 
%                   or 
%              dimVh x 3 matrix \int_{Gamma} sigma n \phi_i
%     MBd    : boundary mass matrix int_{Gamma}lambda \phi_i \phi_j for LHS
%               (only when lambda is specified)
%
% Last modified: November 21, 2017

if ~iscell(g)
    g={g};
end
param=length(g(:));   % 1, 3, 6, or 9 (for 3x3 case)
if size(T.neumann,2)==0
    dimVh=dimFEMspace(T,k);
    switch param
        case {6,9}
            gh=zeros(dimVh,3);
        case {1,3}
            gh=zeros(dimVh,1);
    end
    MBd=sparse(dimVh,dimVh);
    return
end

neuFaces=find(T.faces(4,:)==2);
switch param
    case 1
        testmode='sc';
    case {3,9}
        testmode='fl';
    case 6
        testmode='tr';
end
if nargin==4
    lambda=varargin{1};
else
    lambda=0;
end
if param==9
    [gh(:,1),~,~,MBd]=bdWork3D({g{1,:}},T,k,neuFaces,'test',testmode,lambda);
    gh(:,2)=bdWork3D({g{2,:}},T,k,neuFaces,'test',testmode,lambda);
    gh(:,3)=bdWork3D({g{3,:}},T,k,neuFaces,'test',testmode,lambda);
else
    [gh,~,~,MBd]=bdWork3D(g,T,k,neuFaces,'test',testmode,lambda);
end

return