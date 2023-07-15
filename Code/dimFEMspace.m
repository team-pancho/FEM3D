function d=dimFEMspace(T,k)

% d=dimFEMspace(T,k)
%
% Input: 
%    T    :  tetrahedrization
%    k    :  polynomial degree
% Output:
%    d    :  dimension of the Pk FEM space
% Last modified: July 8, 2016

if ~isfield(T,'edges')
    T=edgesAndFaces(T);
end

Nnod=size(T.coordinates,2);
Nelt=size(T.elements,2);
Nedg=size(T.edges,2);
Nfac=size(T.faces,2);

dimEdges=(k-1)*(k>1);
dimFaces=((k-1)*(k-2))/2*(k>2);
dimInt  =((k-1)*(k-2)*(k-3))/6*(k>3);

d = Nnod+Nedg*dimEdges+Nfac*dimFaces+Nelt*dimInt;

return