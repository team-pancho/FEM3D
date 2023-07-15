function [Tnew,SubT]=subtet(T,Elts)

% [Tnew,SubT]=subtet(T,Elts)
%
% Input:
%     T    : fully enhanced tetrahedrization
%     Elts : a list of integers which represents the 'marked' elements
% Output:
%     Tnew : fully enhanced tetrahedrization
%     SubT : subtet data structure
%             SubT.V : (row) list of marked vertices
%             SubT.E : (row) list of marked edges
%             SubT.F : (row) list of marked faces
%             SubT.EL: (row) list of marked elements
%
% Last modified: October 28, 2016

V = T.elements(:,Elts);
SubT.V = unique(V(:)');
E = abs(T.edgebyelt(:,Elts));
SubT.E = unique(E(:)');
F = T.facebyelt(:,Elts);
SubT.F = unique(F(:)');
SubT.EL = sort(Elts);

Nvert = size(T.coordinates,2);
Nedg  = size(T.edges,2);
Nface = size(T.faces,2);
Nelts = size(T.elements,2);

NNvert = length(SubT.V);
NNedg  = length(SubT.E);
NNface = length(SubT.F);
NNelts = length(SubT.EL);

Tnew.coordinates=T.coordinates(:,SubT.V);
Tnew.elements   =T.elements(:,SubT.EL);
Tnew.edges      =T.edges(:,SubT.E);
Tnew.faces      =T.faces(:,SubT.F);
Tnew.facebyelt  =T.facebyelt(:,SubT.EL);
Tnew.edgebyelt  =abs(T.edgebyelt(:,SubT.EL));
Tnew.faceorient =T.faceorient(:,SubT.EL);
   
transV=zeros(1,Nvert);
transV(SubT.V)=1:NNvert;
transE=zeros(1,Nedg);
transE(SubT.E)=1:NNedg;
transF=zeros(1,Nface);
transF(SubT.F)=1:NNface;

Tnew.elements(:) =transV(Tnew.elements(:));
Tnew.edges(1:2,:)=transV(Tnew.edges(1:2,:));
Tnew.faces(1:3,:)=transV(Tnew.faces(1:3,:));
Tnew.facebyelt(:)=transF(Tnew.facebyelt(:));
Tnew.edgebyelt(:)=transE(Tnew.edgebyelt(:));
Tnew.edgebyelt   =sign(T.edgebyelt(:,SubT.EL)).*Tnew.edgebyelt;

Tnew.dirichlet = Tnew.faces(1:3, find(Tnew.faces(4,:) == 1));
Tnew.neumann   = Tnew.faces(1:3, find(Tnew.faces(4,:) == 2));

Tnew=enhanceGrid3D(Tnew);

return