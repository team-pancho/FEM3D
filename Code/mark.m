function [ SubT ] = mark( T, Elts )

%function [ SubT ] = mark( T, Elts )
%   
% Input:
%       T    : an enhanced tetrahedrization
%       Elts : a list of integers which represents the 'marked' elements
% 
% Output:
%       SubT: a data structure containing the following fields
%             SubT.V : a list of marked vertices
%             SubT.E : a list of marked edges
%             SubT.F : a list of marked faces
%             SubT.EL: a list of marked elements
%             (all lists are sorted to be in increasing order)
% 
% Last Modified: September 9, 2016

V = T.elements(:,Elts);
SubT.V = unique(V(:)');
E = abs(T.edgebyelt(:,Elts));
SubT.E = unique(E(:)');
F = T.facebyelt(:,Elts);
SubT.F = unique(F(:)');
SubT.EL = sort(Elts);

end

