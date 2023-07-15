function eh=errorDFEM3D(u,uh,T,k)

% eh=errorDFEM3D(u,uh,T,k)
% 
% Input:
%    u    : vectorized function of three variables
%    uh   : dim P_k-1 X NElts vector for a discontinuous function
%    T    : basic triangulation
%    k    : polynomial degree of continuous space
%           CAUTION: we really need k-1, but we input k!
% Output:
%    eh   : L2 error
%    
% Last modified: July 8, 2016

T=edgesAndFaces(T);
T=enhanceGrid3D(T);

formula=quadratureFEM(4*k-4,3);
x=T.coordinates(1,:); x=formula(:,1:4)*x(T.elements);
y=T.coordinates(2,:); y=formula(:,1:4)*y(T.elements);
z=T.coordinates(3,:); z=formula(:,1:4)*z(T.elements);
U=u(x,y,z);

Nelts = size(T.elements,2);
dkM1 = size(uh,1)/Nelts;
Uh=reshape(uh, dkM1, Nelts);
uh=bernstein3D(formula(:,2),formula(:,3),formula(:,4),k-1)*Uh;    

eh=sqrt(formula(:,5)'*abs(U-uh).^2*T.volume.');
return
