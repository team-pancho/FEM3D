function sigma=stressAtBaryc(mu,lam,uhtot,T,k)

% sigma=stressAtBaryc(mu,lambda,uhtot,T,k)
% sigma=stressAtBaryc(mu,lambda,{uhtot, Ntime},T,k)
%
% Input:
%    mu,lam: Vectorized function handles of three variables
%    uh   : dim FE-vector for a FE function
%    T    : basic triangulation
%    k    : polynomial degree
% Output:
%    sigma: Nelt x Ntime sized matrix valued stress
%           at barycenter of elements at each time steps
%           
% Last modified: July 1, 2016 by HE

Ne = size(T.elements,2);

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

% Calculating the coordinates of the barycenter
baryc=1/4*(T.coordinates(:,T.elements(1,:))+...
    T.coordinates(:,T.elements(2,:))+...
    T.coordinates(:,T.elements(3,:))+...
    T.coordinates(:,T.elements(4,:)));
barycx=baryc(1,:).'; barycy=baryc(2,:).'; barycz=baryc(3,:).';
[Px,Py,Pz]=bernsteinDer3D(1/4,1/4,1/4,k);

% Decomposing uhtot for several cases
if iscell(uhtot)
    switch length(uhtot)
        case 2
            uhtotal = uhtot{1};
            Ntime = uhtot{2};
            t=2*pi*(1:Ntime)./Ntime;
            len = size(uhtotal,1)/3;
            u = uhtotal(1:len,:);
            v = uhtotal(len+1:2*len,:);
            w = uhtotal(2*len+1:end,:);
            uh=bsxfun(@times,real(u),cos(t))+bsxfun(@times,imag(u),sin(t));
            vh=bsxfun(@times,real(v),cos(t))+bsxfun(@times,imag(v),sin(t));
            wh=bsxfun(@times,real(w),cos(t))+bsxfun(@times,imag(w),sin(t));
        case 3
            uh = uhtot{1};
            vh = uhtot{2};
            wh = uhtot{3};
            Ntime = size(uh,2);
        case 4
            u = uhtot{1};
            v = uhtot{2};
            w = uhtot{3};
            Ntime = uhtot{4};
            t=2*pi*(1:Ntime)./Ntime;
            uh=bsxfun(@times,real(u),cos(t))+bsxfun(@times,imag(u),sin(t));
            vh=bsxfun(@times,real(v),cos(t))+bsxfun(@times,imag(v),sin(t));
            wh=bsxfun(@times,real(w),cos(t))+bsxfun(@times,imag(w),sin(t));
    end
else
    len = size(uhtot,1)/3;
    uh = uhtot(1:len,:);
    vh = uhtot(len+1:2*len,:);
    wh = uhtot(2*len+1:end,:);
    Ntime = size(uh,2);
end
dof = DOF3D(T,k);

% Calculating the partial derivatives at each time step
uhx = zeros(Ne,Ntime);
uhy = zeros(Ne,Ntime);
uhz = zeros(Ne,Ntime);
vhx = zeros(Ne,Ntime);
vhy = zeros(Ne,Ntime);
vhz = zeros(Ne,Ntime);
whx = zeros(Ne,Ntime);
why = zeros(Ne,Ntime);
whz = zeros(Ne,Ntime);

for j=1:Ntime
    Uh=uh(:,j); Uh = Uh(dof);
    Vh=vh(:,j); Vh = Vh(dof);
    Wh=wh(:,j); Wh = Wh(dof);
    
    uhx(:,j)=((Px*bsxfun(@times,c11,Uh)+Py*bsxfun(@times,c21,Uh)+...
        Pz*bsxfun(@times,c31,Uh))).';
    uhy(:,j)=((Px*bsxfun(@times,c12,Uh)+Py*bsxfun(@times,c22,Uh)+...
        Pz*bsxfun(@times,c32,Uh))).';
    uhz(:,j)=((Px*bsxfun(@times,c13,Uh)+Py*bsxfun(@times,c23,Uh)+...
        Pz*bsxfun(@times,c33,Uh))).';
    
    vhx(:,j)=((Px*bsxfun(@times,c11,Vh)+Py*bsxfun(@times,c21,Vh)+...
        Pz*bsxfun(@times,c31,Vh))).';
    vhy(:,j)=((Px*bsxfun(@times,c12,Vh)+Py*bsxfun(@times,c22,Vh)+...
        Pz*bsxfun(@times,c32,Vh))).';
    vhz(:,j)=((Px*bsxfun(@times,c13,Vh)+Py*bsxfun(@times,c23,Vh)+...
        Pz*bsxfun(@times,c33,Vh))).';
    
    whx(:,j)=((Px*bsxfun(@times,c11,Wh)+Py*bsxfun(@times,c21,Wh)+...
        Pz*bsxfun(@times,c31,Wh))).';
    why(:,j)=((Px*bsxfun(@times,c12,Wh)+Py*bsxfun(@times,c22,Wh)+...
        Pz*bsxfun(@times,c32,Wh))).';
    whz(:,j)=((Px*bsxfun(@times,c13,Wh)+Py*bsxfun(@times,c23,Wh)+...
        Pz*bsxfun(@times,c33,Wh))).';
end

% Using the partial derivatives calculating the stress
mubaryc = mu(barycx,barycy,barycz);
lambdabaryc = lam(barycx,barycy,barycz);
div=uhx+vhy+whz;

sxx=2*bsxfun(@times,mubaryc,uhx) + bsxfun(@times,lambdabaryc,div);
syy=2*bsxfun(@times,mubaryc,vhy) + bsxfun(@times,lambdabaryc,div);
szz=2*bsxfun(@times,mubaryc,whz) + bsxfun(@times,lambdabaryc,div);

sxy=bsxfun(@times,mubaryc,(vhx+uhy));
sxz=bsxfun(@times,mubaryc,(whx+uhz));
syz=bsxfun(@times,mubaryc,(why+vhz));


sigma=sqrt(sxx.^2+syy.^2+szz.^2+2*sxy.^2+2*sxz.^2+2*syz.^2);

return
