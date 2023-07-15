function [sigmaxx,sigmayy,sigmazz,sigmayz,sigmaxz,sigmaxy]=...
                      stressPostprocessing(uh,vh,wh,lambda,mu,T,k,varargin)
% [sigmaxx,sigmayy,sigmazz,sigmayz,sigmaxz,sigmaxy]=...
%                     stressPostprocessing(uh,vh,wh,lambda,mu,T,k)
% [sigmaxx,sigmayy,sigmazz,sigmayz,sigmaxz,sigmaxy]=...
%                     stressPostprocessing(uh,vh,wh,lambda,mu,T,k,tag)
%
% Input:
%    uh          : dim FE-vector for a FE function (1st component)
%    vh          : dim FE-vector for a FE function (2nd component)
%    wh          : dim FE-vector for a FE function (3rd component)
%    lambda      : vectorized function of three variables
%    mu          : vectorized function of three variables
%    T           : basic triangulation
%    k           : polynomial degree
%    tag         : an indicator that you want the averaged stress on each
%                  element rather than the value on nodes (can be anything)
%
% Output:
%    sigmaxx,... : dimP_{k-1}*Nelt x 1. Each represents for coefficients
%                  of stress in piecewise P_{k-1} polynomial space
%                  or 
%                  1 x Nelt. Each is a row vector of 
%                  averaged stress on each element
%
% Last modified: Nov 21, 2017

T=edgesAndFaces(T);
T=enhanceGrid3D(T);
Nelt=size(T.elements,2);
    
formula=quadratureFEM(4*k,3);
x=T.coordinates(1,:); x=formula(:,1:4)*x(T.elements);
y=T.coordinates(2,:); y=formula(:,1:4)*y(T.elements);
z=T.coordinates(3,:); z=formula(:,1:4)*z(T.elements);

Mu=mu(x,y,z);
La=lambda(x,y,z);

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

% Evaluation of all stresses at all quadrature points

Uh=uh(DOF3D(T,k));
Vh=vh(DOF3D(T,k));
Wh=wh(DOF3D(T,k));
[Px,Py,Pz]=bernsteinDer3D(formula(:,2),formula(:,3),formula(:,4),k);

uhx=Px*bsxfun(@times,c11,Uh)+Py*bsxfun(@times,c21,Uh)...
        +Pz*bsxfun(@times,c31,Uh);
uhy=Px*bsxfun(@times,c12,Uh)+Py*bsxfun(@times,c22,Uh)...
        +Pz*bsxfun(@times,c32,Uh);
uhz=Px*bsxfun(@times,c13,Uh)+Py*bsxfun(@times,c23,Uh)...
        +Pz*bsxfun(@times,c33,Uh);

vhx=Px*bsxfun(@times,c11,Vh)+Py*bsxfun(@times,c21,Vh)...
        +Pz*bsxfun(@times,c31,Vh);
vhy=Px*bsxfun(@times,c12,Vh)+Py*bsxfun(@times,c22,Vh)...
        +Pz*bsxfun(@times,c32,Vh);
vhz=Px*bsxfun(@times,c13,Vh)+Py*bsxfun(@times,c23,Vh)...
        +Pz*bsxfun(@times,c33,Vh);

whx=Px*bsxfun(@times,c11,Wh)+Py*bsxfun(@times,c21,Wh)...
        +Pz*bsxfun(@times,c31,Wh);
why=Px*bsxfun(@times,c12,Wh)+Py*bsxfun(@times,c22,Wh)...
        +Pz*bsxfun(@times,c32,Wh);
whz=Px*bsxfun(@times,c13,Wh)+Py*bsxfun(@times,c23,Wh)...
        +Pz*bsxfun(@times,c33,Wh);

sigmaxx=2*Mu.*uhx+La.*(uhx+vhy+whz);
sigmayy=2*Mu.*vhy+La.*(uhx+vhy+whz);
sigmazz=2*Mu.*whz+La.*(uhx+vhy+whz);
sigmaxy=Mu.*(uhy+vhx);
sigmaxz=Mu.*(uhz+whx);
sigmayz=Mu.*(vhz+why);

Sigmas=[sigmaxx,sigmayy,sigmazz,sigmaxy,sigmaxz,sigmayz];

% L2 projection element by element

R =bernstein3D(formula(:,2),formula(:,3),formula(:,4),k-1);
wR=bsxfun(@times,formula(:,5),R);
SigmaProj=(R'*wR)\(wR'*Sigmas);

if nargin > 7
    SigmaProj = sum(wR*SigmaProj);
    sigmaxx = SigmaProj(1:Nelt);
    sigmayy=SigmaProj(Nelt+1:2*Nelt);
    sigmazz=SigmaProj(2*Nelt+1:3*Nelt);
    sigmaxy=SigmaProj(3*Nelt+1:4*Nelt);
    sigmaxz=SigmaProj(4*Nelt+1:5*Nelt);
    sigmayz=SigmaProj(5*Nelt+1:6*Nelt);
else
    sigmaxx=SigmaProj(:,1:Nelt);
    sigmayy=SigmaProj(:,Nelt+1:2*Nelt);
    sigmazz=SigmaProj(:,2*Nelt+1:3*Nelt);
    sigmaxy=SigmaProj(:,3*Nelt+1:4*Nelt);
    sigmaxz=SigmaProj(:,4*Nelt+1:5*Nelt);
    sigmayz=SigmaProj(:,5*Nelt+1:6*Nelt);
    sigmaxx=sigmaxx(:);
    sigmayy=sigmayy(:);
    sigmazz=sigmazz(:);
    sigmaxy=sigmaxy(:);
    sigmaxz=sigmaxz(:);
    sigmayz=sigmayz(:);
end


return
