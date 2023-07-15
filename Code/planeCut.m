function Elts = planeCut( T, planes )

% Elts = planeCut( T, planes )   
% 
% Input:
%       T       : an enhanced tetrahedrization
%       planes  : a number of planes X 4 where each row contains the
%                 coefficients of a plane 
% 
% Output:
%       Elts    : the 'marked' elements whose barycenters are above all of
%                 the planes
% 
% Last Modified: September 13, 2016

bary = (T.coordinates(:, T.elements(1,:))...
          + T.coordinates(:, T.elements(2,:))...
          + T.coordinates(:, T.elements(3,:))...
          + T.coordinates(:, T.elements(4,:)))/4;

bary = [bary; ones(1, size(bary,2))];

test = planes*bary;
test = test>0;
Elts = prod(test,1);
Elts = find(Elts==1);

end

