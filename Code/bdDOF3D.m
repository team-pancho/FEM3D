function [DirDOF,NeuDOF,NeuAssem,dirFaces,neuFaces,dirlist,free] = bdDOF3D(T,k)

% [DirDOF,NeuDOF,NeuAssem,dirFaces,neuFaces,dirlist,free] = bdDOF3D(T,k)
%
% Input:
% T             : Data structure. Enhanced 3D FEM tetrahedrization
% k             : Polynomial degree
% Output:
% DirDOF        : dim P_k(F) x NDirFaces. Matrix
%                 containing global boundary DOF per Dirichlet face counted 
%                 columnwise.
% NeuDOF        : dim P_k(F) x NNeuFaces. Matrix
%                 containing global boundary DOF per Neumann face counted 
%                 columnwise.
% NeuAssem      : dim P_k(F) x dim P_k(F) x NneuFaces. Neumann-assembly matrix
% dirFaces      : 1 x NdirFaces. Vector containing the global numbering of
%                 the Dirichlet faces.
% neuFaces      : 1 x NneuFaces. Vector containing the global numbering of
%                 the Neumann faces.
% dirlist       : NdirDOF x 1. List of Dirichlet DOF
% free          : Nfree x 1. vector containing the global numbering of the
%                 free DOF.
%
% NOTE: DOF on Dirichlet/Neumann transition nodes or edges are counted
% twice, appearing in both DirDOF and NeuDOF.
%
% Last modified: January 18, 2017

dirFaces = find(T.faces(4,:)==1);
neuFaces = find(T.faces(4,:)==2);
[DirDOF,~] = computeBDDOF3D(T,k,dirFaces);
[NeuDOF,NeuAssem] = computeBDDOF3D(T,k,neuFaces);

dirlist=unique(DirDOF(:));
free   =(1:dimFEMspace(T,k))'; free(dirlist)=[];

return



