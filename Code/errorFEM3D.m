function [eh,hh]=errorFEM3D(u,uh,T,k)

% eh=errorFEM3D(u,uh,T,k)
% [eh,hh]=errorFEM3D({u,ux,uy,uz},uh,T,k)
%
% Input:
%    u    : vectorized function of two variables, or
%    u,ux,uy,uz  : vectorized functions of three variables
%    uh   : dim FE-vector for a FE function
%    T    : basic triangulation
%    k    : polynomial degree
% Output:
%    eh   : L2 error
%    hh   : H1 error
% Last modified: May 1, 2015

T=edgesAndFaces(T);
T=enhanceGrid3D(T);
formula=quadratureFEM(4*k,3);
x=T.coordinates(1,:); x=formula(:,1:4)*x(T.elements);
y=T.coordinates(2,:); y=formula(:,1:4)*y(T.elements);
z=T.coordinates(3,:); z=formula(:,1:4)*z(T.elements);
if iscell(u)
    U=u{1}(x,y,z);
    Ux=u{2}(x,y,z);
    Uy=u{3}(x,y,z);
    Uz=u{4}(x,y,z);
else
    U=u(x,y,z);
end
Uh=uh(DOF3D(T,k));
uh=bernstein3D(formula(:,2),formula(:,3),formula(:,4),k)*Uh;    
eh=sqrt(formula(:,5)'*abs(U-uh).^2*T.volume.');
if ~iscell(u)
    hh=0;
    return
end

% Geometric coefficients for change of variables

x12=T.coordinates(1,T.elements(2,:))-T.coordinates(1,T.elements(1,:));  
y12=T.coordinates(2,T.elements(2,:))-T.coordinates(2,T.elements(1,:));  
z12=T.coordinates(3,T.elements(2,:))-T.coordinates(3,T.elements(1,:)); 
x13=T.coordinates(1,T.elements(3,:))-T.coordinates(1,T.elements(1,:));   
y13=T.coordinates(2,T.elements(3,:))-T.coordinates(2,T.elements(1,:));
z13=T.coordinates(3,T.elements(3,:))-T.coordinates(3,T.elements(1,:));
x14=T.coordinates(1,T.elements(4,:))-T.coordinates(1,T.elements(1,:));   
y14=T.coordinates(2,T.elements(4,:))-T.coordinates(2,T.elements(1,:));
z14=T.coordinates(3,T.elements(4,:))-T.coordinates(3,T.elements(1,:));

% Entries of CK

sqdet=1./(6*T.volume);
c11=sqdet.*(y13.*z14 - y14.*z13);
c12=sqdet.*(x14.*z13 - x13.*z14);
c13=sqdet.*(x13.*y14 - x14.*y13);
c21=sqdet.*(y14.*z12 - y12.*z14);
c22=sqdet.*(x12.*z14 - x14.*z12);
c23=sqdet.*(x14.*y12 - x12.*y14);
c31=sqdet.*(y12.*z13 - z12.*y13);
c32=sqdet.*(x13.*z12 - x12.*z13);
c33=sqdet.*(x12.*y13 - x13.*y12);


[Px,Py,Pz]=bernsteinDer3D(formula(:,2),formula(:,3),formula(:,4),k);
uhx=Px*bsxfun(@times,c11,Uh)+Py*bsxfun(@times,c21,Uh)+Pz*bsxfun(@times,c31,Uh);
uhy=Px*bsxfun(@times,c12,Uh)+Py*bsxfun(@times,c22,Uh)+Pz*bsxfun(@times,c32,Uh);
uhz=Px*bsxfun(@times,c13,Uh)+Py*bsxfun(@times,c23,Uh)+Pz*bsxfun(@times,c33,Uh);
hh=sqrt(formula(:,5)'*(abs(Ux-uhx).^2+abs(Uy-uhy).^2 +abs(Uz-uhz).^2)*T.volume.');

return
