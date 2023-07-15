function [ t, wts ] = gaussJacobi(nt,alpha)

%function [ t, wts ] = cdgqf( nt,alpha)
%   
%  Adapted from code written by John Burkardt, which in turn was adapted 
%  from code written by Sylvan Elhay and Jaroslav Kautsky. 
% 
% Inputs:
%       nt    : number of desired points
%       alpha : power of (1-t) as in \int_{-1}^1 f(t) (1-t)^alpha dt
%
% Outputs:
%       t     : quadrature nodes
%       wts   : quadrature weights  (both for Gauss-Jacobi)
%
% Last Modified: March 24, 2016
  
bj = zeros(nt,1);
aj = zeros(nt,1);
  
aj(1) = -alpha / (alpha +2);
bj(1) = sqrt(4*(1 + alpha)/((alpha + 3)*(alpha + 2)^2));

aj(2:nt) = -alpha^2./((2*(2:nt) + alpha - 2).*(2*(2:nt) + alpha));
bj(2:nt) = sqrt((2*(2:nt).*((2:nt) + alpha)).^2 ...
           ./(((2*(2:nt) + alpha).^2 -1).*(2*(2:nt) + alpha).^2));

%   Diagonalize the Jacobi matrix.
  
A = diag(aj) + diag(bj(1:nt-1),-1) + diag(bj(1:nt-1),1);
[V,D] = eig(A);

%  Compute the knots and weights.

t = diag(D);

wts = zeros(nt, 1);
wts(1) = sqrt(2^(alpha + 1)*gamma(alpha + 1)*1/gamma(alpha +2));
wts = V'*wts;
wts = wts.^2;

  return
end
  