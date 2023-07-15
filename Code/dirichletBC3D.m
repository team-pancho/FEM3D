function [uh,dir,free]=dirichletBC3D(u,T,k)

% [uh,dir,free]=dirichletBC3D(u,T,k)
%
% Input:
%     u    : vectorized function of three variables 
%              or cell array with M vectorized functions of three vars
%     T    : enhanced triangulation
%     k    : polynomial degree
% Output:
%    uh    : dim FE_h column vector with assigned Dirichlet DOF
%               or dim FE_h x M matrix
%    dir   : List of Dirichlet degrees of freedom
%    free  : List of ALL non Dirichlet DOF
%
% Last modified: January 24, 2017

dimVh=dimFEMspace(T,k);
if size(T.dirichlet,2)==0    
    if iscell(u)
        uh=zeros(dimVh,length(u));
    else
    uh=zeros(dimVh,1);
    end
    dir=unique([]); free=(1:dimVh)';
    return
end
dirFaces=find(T.faces(4,:)==1);
if isa(u,'function_handle')
    u={u};
end
[uh,~,dirDOF,~]=bdWork3D(u,T,k,dirFaces,'proj','sc',0);
dir  =unique(dirDOF(:));
free =(1:dimVh)'; free(dir)=[];

return