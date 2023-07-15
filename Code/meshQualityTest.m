function [ dF,dL,dV ] = meshQualityTest( T )

%function [ dF,dL,dV ] = meshQualityTest( T )
% 
% Inputs:
%       T : a basic tetrahedrization
% 
% Outputs:
%       dF : maximum face deformation  
%       dL : maximum edge deformation
%       dV : ratio of largest/smallest volumes of elements
% 
% Last Modified: February 19, 2016

T = edgesAndFaces(T);
T = enhanceGrid3D(T);

F = T.area(T.facebyelt);

lengths = sqrt(sum((T.coordinates(:, T.edges(1,:))...
                     - T.coordinates(:, T.edges(2,:))).^2));

L = lengths(abs(T.edgebyelt));

dF = max(max(F)/min(F));
dL = max(max(L)/min(L));
dV = max(T.volume)/min(T.volume);

end

