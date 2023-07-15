function Where=embed(Tbig,SubT,k)

% Where=embed(Tbig,SubT,k)
%
% Input:
%    Tbig   : full enhanced tetrahedrization
%    SubT   : SubTet data structure 
%    k      : polynomial degree
% Output:
%    Where  : dimVh(Tsmall) column vector with embedding operator from
%               Vh(Tsmall) to Vh(Tbig)
% 
% Last modified: September 16, 2016

Nvert=size(Tbig.coordinates,2);
Nedge=size(Tbig.edges,2);
Nface=size(Tbig.faces,2);

dimE  =k-1;
dimF  =(k-1)*(k-2)/2;
dimEL =(k-1)*(k-2)*(k-3)/6;

Vert =SubT.V;
Edges=bsxfun(@plus,(1:dimE)',dimE*(SubT.E-1));
Faces=bsxfun(@plus,(1:dimF)',dimF*(SubT.F-1));
Elts =bsxfun(@plus,(1:dimEL)',dimEL*(SubT.EL-1));
Where=[Vert(:);...
       Edges(:)+Nvert;...
       Faces(:)+Nvert+Nedge*dimE;...
       Elts(:)+Nvert+Nedge*dimE+Nface*dimF];
return