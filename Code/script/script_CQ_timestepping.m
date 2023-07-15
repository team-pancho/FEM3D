% A script to figure out CQ time stepping with scalar problems with 3D FEM.
% Last modified: May 8, 2016

if ~exist('external')
    k=input('polynomial degree: ');
    TT=input('final time: ');
    M=input('number of time steps: ');
    lev = input('refinement level: ');
end

kappa=TT/M;
DBC=[1 2 3 4 5 6]; %dirichlet faces

T=meshCube(lev,lev,lev,DBC);
T=edgesAndFaces(T);
T=enhanceGrid3D(T);

% interior planewave and derivatives
HH =  @(x) x.^5.*(1-5*(x-1)+15*(x-1).^2-35*(x-1).^3+70*(x-1).^4-...
           126*(x-1).^5).*(x>0).*(x<1)+(x>=1);
HHp = @(x) (5*x.^4.*...
           (1-5*(x-1)+15*(x-1).^2-35*(x-1).^3+70*(x-1).^4-126*(x-1).^5)...
           +x.^5.*(-5+30*(x-1)-105*(x-1).^2+280*(x-1).^3-630*(x-1).^4))...
           .*(x>0).*(x<1);
HHpp= @(x) -1260*(x-1).^4.*(9*x-4).*(x.^3).*(x>0).*(x<1);          

alpha=1; ll=3;   % [0,alpha] - [alpha,l-alpha] [l-alpha l]
window  = @(x) HH(x/alpha).*HH((ll-x)./alpha);
windowp = @(x) 1/alpha*(HHp(x/alpha).*HH((ll-x)./alpha)...
                        -HH(x/alpha).*HHp((ll-x)./alpha));
windowpp= @(x) 1/alpha^2*(HH((ll-x)./alpha).*HHpp(x/alpha)+...
                HH(x/alpha).*HHpp((ll-x)./alpha) ...
                -2*HHp((ll-x)./alpha).*HH(x/alpha));

signal=@(t) sin(2*t).*window(t);
signalp=@(t) sin(2*t).*windowp(t)+2*cos(2*t).*window(t);
signalpp=@(t) sin(2*t).*windowpp(t)+4*cos(2*t).*windowp(t)-4*sin(2*t).*window(t);

d=[1 1 1]/sqrt(3);
tlag=1;

u=@(x1,x2,x3,t) signal(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
ux=@(x1,x2,x3,t) -d(1)*signalp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uy=@(x1,x2,x3,t) -d(2)*signalp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uz=@(x1,x2,x3,t) -d(3)*signalp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uxx=@(x1,x2,x3,t) d(1)^2*signalpp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uyy=@(x1,x2,x3,t) d(2)^2*signalpp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uzz=@(x1,x2,x3,t) d(3)^2*signalpp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uxy=@(x1,x2,x3,t) d(1)*d(2)*signalpp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uxz=@(x1,x2,x3,t) d(1)*d(3)*signalpp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uyz=@(x1,x2,x3,t) d(2)*d(3)*signalpp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);
uyx=@(x1,x2,x3,t) uxy(x1,x2,x3,t);
uzx=@(x1,x2,x3,t) uxz(x1,x2,x3,t);
uzy=@(x1,x2,x3,t) uyz(x1,x2,x3,t);
utt=@(x1,x2,x3,t) signalpp(-(x1*d(1)+x2*d(2)+x3*d(3))+t-tlag);

% FEM parameters (wavespeed and material properties)
c=@(x,y,z) 1+0*x;
kxx = @(x,y,z) 1+0.5*(x.^2+y.^2+z.^2);
kxy = @(x,y,z) 0.25+0.5*(x.^2+y.^2+z.^2);
kxz = @(x,y,z) 2+0.5*(x.^2+y.^2+z.^2);
kyy = @(x,y,z) 3+0.5*(x.^2+y.^2+z.^2);
kyz = @(x,y,z) 0.5+0.5*(x.^2+y.^2+z.^2);
kzz = @(x,y,z) 6+0.5*(x.^2+y.^2+z.^2);
dxkxx = @(x,y,z) x;
dxkxy = @(x,y,z) x; 
dxkxz = @(x,y,z) x;
dykxy = @(x,y,z) y;
dykyy = @(x,y,z) y;
dykyz = @(x,y,z) y;
dzkzz = @(x,y,z) z;
dzkxz = @(x,y,z) z;
dzkyz = @(x,y,z) z;

% matrices 
Mh = massMatrix3D(c,T,k);
Sh=stiffnessMatrices3D(kxx,kxy,kxz,kyy,kyz,kzz,T,k);
Sh=Sh{1,1}+Sh{2,2}+Sh{3,3}+Sh{1,2}+Sh{1,2}'+Sh{1,3}+Sh{1,3}'+Sh{2,3}+Sh{2,3}';

dimVh = size(Mh,1);

% \kappa \grad u 
qx = @(x,y,z,t) kxx(x,y,z).*ux(x,y,z,t)+kxy(x,y,z).*uy(x,y,z,t)+kxz(x,y,z).*uz(x,y,z,t);
qy = @(x,y,z,t) kxy(x,y,z).*ux(x,y,z,t)+kyy(x,y,z).*uy(x,y,z,t)+kyz(x,y,z).*uz(x,y,z,t);
qz = @(x,y,z,t) kxz(x,y,z).*ux(x,y,z,t)+kyz(x,y,z).*uy(x,y,z,t)+kzz(x,y,z).*uz(x,y,z,t);

% RHS such that the solution satisfies the wave equation
f=@(x,y,z,t) -(kxx(x,y,z).*uxx(x,y,z,t)+dxkxx(x,y,z).*ux(x,y,z,t)+...
             kxy(x,y,z).*uxy(x,y,z,t)+dxkxy(x,y,z).*uy(x,y,z,t)+...
             kxz(x,y,z).*uxz(x,y,z,t)+dxkxz(x,y,z).*uz(x,y,z,t))...
             -(kxy(x,y,z).*uxy(x,y,z,t)+dykxy(x,y,z).*ux(x,y,z,t)+...
             kyy(x,y,z).*uyy(x,y,z,t)+dykyy(x,y,z).*uy(x,y,z,t)+...
             kyz(x,y,z).*uyz(x,y,z,t)+dykyz(x,y,z).*uz(x,y,z,t))...
             -(kxz(x,y,z).*uxz(x,y,z,t)+dzkxz(x,y,z).*ux(x,y,z,t)+...
             kyz(x,y,z).*uyz(x,y,z,t)+dzkyz(x,y,z).*uy(x,y,z,t)+...
             kzz(x,y,z).*uzz(x,y,z,t)+dzkzz(x,y,z).*uz(x,y,z,t))...
             +utt(x,y,z,t);

neu = zeros(dimVh,M+1);
uh  = zeros(dimVh,M+1);
rhs = zeros(dimVh,M+1);

[~,~,~,~,~,dir,free] = bdDOF3D(T,k);

%the parfor helps a good deal
parfor t=0:M 
    wx=@(x,y,z) qx(x,y,z,t*kappa);
    wy=@(x,y,z) qy(x,y,z,t*kappa);
    wz=@(x,y,z) qz(x,y,z,t*kappa);
    neu(:,t+1) = neumannBC3D({wx,wy,wz},T,k);
    uh(:,t+1)  = dirichletBC3D(@(x,y,z) u(x,y,z,t*kappa),T,k);
    rhs(:,t+1) = loadVector3D(@(x,y,z) f(x,y,z,t*kappa),T,k);
end

RHS = rhs+neu;

FEMsolve = @(s,v,rhs) (s^2*Mh(free,free)+Sh(free,free))\...
    (rhs(free)-s^2*Mh(free,dir)*v(dir)-Sh(free,dir)*v(dir));

% CQ stuff
p=@(z) 2*(1-z)/(1+z);

omega = exp(2*pi*1i/(M+1));
R = eps^(0.5/(M+1));

uh = bsxfun(@times,uh,R.^(0:M));    
uh = fft(uh,[],2);               % DFT by columns (\hat H)

h = bsxfun(@times,RHS,R.^(0:M));    
h = fft(h,[],2);               % DFT by columns (\hat H)

% not yet parallelizable
for l=0:floor((M+1)/2)
    s=p(R*omega^(-l))/kappa;
    uh(free,l+1)=FEMsolve(s,uh(:,l+1),h(:,l+1));    % \hat v   
end
uh(:,M+2-(1:floor(M/2)))=conj(uh(:,2:floor(M/2)+1));

uh=real(ifft(uh,[],2));                           % v
uh=bsxfun(@times,uh,R.^(-(0:M)));

% exact solution at final time
u =@(x,y,z)  u(x,y,z,TT);
ux=@(x,y,z) ux(x,y,z,TT);
uy=@(x,y,z) uy(x,y,z,TT);
uz=@(x,y,z) uz(x,y,z,TT);

% computation of errors

[eL2,eH1]=errorFEM3D({u,ux,uy,uz},uh(:,end),T,k);

if ~exist('external')
    disp(eL2);
    disp(eH1);
else
    errors = [errors; eL2 eH1];
end


