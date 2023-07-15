function [listDOF,assembleDOF]=computeBDDOF3D(T,k,list)

% [listDOF,asseDOF]=computeBDDOF3D(T,k,list)
%  
% Input:
%    T       : full tetrahedrization
%    k       : polynomial degree
%    list    : row vector with numbers from 1 to Nfac
% Output:
%    listDOF : dim P_k(F) x #list matrix with global DOF for listed faces
%    asseDOF : dim P_k(F) x dim P_k(F) x #list COLUMN assembly array
%
% Last modified: January 18, 2017

Nnod = size(T.coordinates,2);
Nedg = size(T.edges,2);

%  Edgebyface 
Edges = T.edges(1:2,:)';
[E1,I1] = sort(T.faces([1 2],:)',2);
[E2,I2] = sort(T.faces([2 3],:)',2);
[E3,I3] = sort(T.faces([3 1],:)',2);
[~,EBF1] = ismember(E1,Edges,'rows');
[~,EBF2] = ismember(E2,Edges,'rows');
[~,EBF3] = ismember(E3,Edges,'rows');
EdgebyFace = [EBF1 EBF2 EBF3]';
Orient = [I1(:,2)-I1(:,1) I2(:,2)-I2(:,1) I3(:,2)-I3(:,1)]';
listEdges = EdgebyFace(:,list); listEdges = listEdges(:)';
listOrient = Orient(:,list); listOrient = listOrient(:)';
    
% Nodal DOF
listnodeDOF=T.faces(1:3,list);

% Edge DOF (orientation enforced in second line)
listedgeDOF = bsxfun(@plus,(k-1)*(listEdges-1) +Nnod,(1:k-1)');
listedgeDOF(:,listOrient<0) = flipud(listedgeDOF(:,listOrient<0)); 
listedgeDOF = reshape(listedgeDOF,3*(k-1),length(list));

% Face DOF (offset by nodal and edge DOF)
dimk=(k-1)*(k-2)/2; % number of internal/face DOF
listfaceDOF = bsxfun(@plus,dimk*(list-1)+Nnod+(k-1)*Nedg,(1:dimk)');
listfaceDOF = reshape(listfaceDOF,dimk,length(list));

listDOF     = [listnodeDOF; listedgeDOF; listfaceDOF];

% Assembly matrix (columns only)

dk = 3+3*(k-1)+(k-1)*(k-2)/2;   % dim P_k(F)
assembleDOF = repmat(listDOF(:)',[dk,1]);
assembleDOF = reshape(assembleDOF,[dk,dk,length(list)]);

end

