% A script to couple CQ time stepping with 3D FEM for elasticity
% Last modified: May 14, 2016

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
e=d;
c=1;
tlag=1;

% material parameters & exact solution
mu  =@(x,y,z) 5+(x.^2+y.^2+z.^2);
lam =@(x,y,z) 2+(x.^2+y.^2+z.^2);
rho =@(x,y,z) 3+(x.^2+y.^2+z.^2);

lamx =@(x,y,z) 2*x;
lamy =@(x,y,z) 2*y;
lamz =@(x,y,z) 2*z;

mux =@(x,y,z) 2*x;
muy =@(x,y,z) 2*y;
muz =@(x,y,z) 2*z;

[u,v,w,ut,vt,wt,utt,vtt,wtt,...
                ux,uy,uz,vx,vy,vz,wx,wy,wz,...
                uxx,uxy,uxz,uyy,uyz,uzz,...
                vxx,vxy,vxz,vyy,vyz,vzz,...
                wxx,wxy,wxz,wyy,wyz,wzz] ...
                = vectorPlaneWave(e,d,c,tlag,{signal,signalp,signalpp});

% rhs so the equation is satisfied
div =@(x,y,z,t) ux(x,y,z,t)+vy(x,y,z,t)+wz(x,y,z,t);
divx=@(x,y,z,t) uxx(x,y,z,t)+vxy(x,y,z,t)+wxz(x,y,z,t);
divy=@(x,y,z,t) uxy(x,y,z,t)+vyy(x,y,z,t)+wyz(x,y,z,t);
divz=@(x,y,z,t) uxz(x,y,z,t)+vyz(x,y,z,t)+wzz(x,y,z,t);

sxx=@(x,y,z,t) 2*mu(x,y,z).*ux(x,y,z,t)+lam(x,y,z).*div(x,y,z,t);
syy=@(x,y,z,t) 2*mu(x,y,z).*vy(x,y,z,t)+lam(x,y,z).*div(x,y,z,t);
szz=@(x,y,z,t) 2*mu(x,y,z).*wz(x,y,z,t)+lam(x,y,z).*div(x,y,z,t);

sxy=@(x,y,z,t) mu(x,y,z).*(vx(x,y,z,t)+uy(x,y,z,t));
sxz=@(x,y,z,t) mu(x,y,z).*(wx(x,y,z,t)+uz(x,y,z,t));
syz=@(x,y,z,t) mu(x,y,z).*(wy(x,y,z,t)+vz(x,y,z,t));

fx=@(x,y,z,t) -(lam(x,y,z).*divx(x,y,z,t)+lamx(x,y,z).*div(x,y,z,t)...
    +2*mu(x,y,z).*uxx(x,y,z,t)+2*mux(x,y,z).*ux(x,y,z,t)...
    +mu(x,y,z).*(vxy(x,y,z,t)+uyy(x,y,z,t))+muy(x,y,z).*(vx(x,y,z,t)+uy(x,y,z,t))...
    +mu(x,y,z).*(wxz(x,y,z,t)+uzz(x,y,z,t))+muz(x,y,z).*(wx(x,y,z,t)+uz(x,y,z,t)))...
    +rho(x,y,z).*utt(x,y,z,t);

fy=@(x,y,z,t) -(lam(x,y,z).*divy(x,y,z,t)+lamy(x,y,z).*div(x,y,z,t)...
    +2*mu(x,y,z).*vyy(x,y,z,t)+2*muy(x,y,z).*vy(x,y,z,t)...
    +mu(x,y,z).*(wyz(x,y,z,t)+vzz(x,y,z,t))+muz(x,y,z).*(wy(x,y,z,t)+vz(x,y,z,t))...
    +mu(x,y,z).*(uxy(x,y,z,t)+vxx(x,y,z,t))+mux(x,y,z).*(uy(x,y,z,t)+vx(x,y,z,t)))...
    +rho(x,y,z).*vtt(x,y,z,t);

fz=@(x,y,z,t) -(lam(x,y,z).*divz(x,y,z,t)+lamz(x,y,z).*div(x,y,z,t)...
    +2*mu(x,y,z).*wzz(x,y,z,t)+2*muz(x,y,z).*wz(x,y,z,t)...
    +mu(x,y,z).*(uxz(x,y,z,t)+wxx(x,y,z,t))+mux(x,y,z).*(uz(x,y,z,t)+wx(x,y,z,t))...
    +mu(x,y,z).*(vyz(x,y,z,t)+wyy(x,y,z,t))+muy(x,y,z).*(vz(x,y,z,t)+wy(x,y,z,t)))...
    +rho(x,y,z).*wtt(x,y,z,t);

% matrices
Sm = stiffnessMatrices3D(mu,mu,mu,mu,mu,mu,T,k);
Sl = stiffnessMatrices3D(lam,lam,lam,lam,lam,lam,T,k);
Mh = massMatrix3D(rho,T,k);

Sh=[2*Sm{1,1}+Sl{1,1}+Sm{2,2}+Sm{3,3}     Sl{1,2}+Sm{1,2}'      Sl{1,3}+Sm{1,3}';...
     Sm{1,2}+Sl{1,2}'      2*Sm{2,2}+Sm{1,1}+Sl{2,2}+Sm{3,3}     Sl{2,3}+Sm{2,3}';...
     Sl{1,3}'+Sm{1,3}      Sl{2,3}'+Sm{2,3}     2*Sm{3,3}+Sm{1,1}+Sm{2,2}+Sl{3,3}];

O=sparse(size(Mh,1),size(Mh,2));
Mh = [Mh O O;...
      O Mh O;...
      O O Mh];   

Ndof  = size(Sm{1},1);
dimVh = 3*Ndof;
uh = zeros(dimVh,M+1);
RHS = zeros(dimVh,M+1);

parfor t=0:M
    uk=@(x,y,z) u(x,y,z,t*kappa);
    vk=@(x,y,z) v(x,y,z,t*kappa);
    wk=@(x,y,z) w(x,y,z,t*kappa);
    uh(:,t+1) = reshape(dirichletBC3D({uk,vk,wk},T,k),dimVh,1);
    fxk=@(x,y,z) fx(x,y,z,t*kappa);
    fyk=@(x,y,z) fy(x,y,z,t*kappa);
    fzk=@(x,y,z) fz(x,y,z,t*kappa);
    RHS(:,t+1) = reshape(loadVector3D({fxk,fyk,fzk},T,k),dimVh,1);
    sigma={};
    sigma{1} = @(x,y,z) sxx(x,y,z,t*kappa);
    sigma{2} = @(x,y,z) sxy(x,y,z,t*kappa);
    sigma{3} = @(x,y,z) sxz(x,y,z,t*kappa);
    sigma{4} = @(x,y,z) syy(x,y,z,t*kappa);
    sigma{5} = @(x,y,z) syz(x,y,z,t*kappa);
    sigma{6} = @(x,y,z) szz(x,y,z,t*kappa);
    RHS(:,t+1) = RHS(:,t+1) + reshape(neumannBC3D(sigma,T,k),dimVh,1);
end

[~,~,~,~,~,dir,free]=bdDOF3D(T,k);

dir  = [dir(:) ; Ndof+dir(:) ;2*Ndof+dir(:)];
free = [free(:); Ndof+free(:);2*Ndof+free(:)];

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

u =@(x,y,z)  u(x,y,z,TT);
ux=@(x,y,z) ux(x,y,z,TT);
uy=@(x,y,z) uy(x,y,z,TT);
uz=@(x,y,z) uz(x,y,z,TT);

v =@(x,y,z)  v(x,y,z,TT);
vx=@(x,y,z) vx(x,y,z,TT);
vy=@(x,y,z) vy(x,y,z,TT);
vz=@(x,y,z) vz(x,y,z,TT);

w =@(x,y,z)  w(x,y,z,TT);
wx=@(x,y,z) wx(x,y,z,TT);
wy=@(x,y,z) wy(x,y,z,TT);
wz=@(x,y,z) wz(x,y,z,TT);

[eL2u,eH1u]=errorFEM3D({u,ux,uy,uz},uh(1:Ndof,end),T,k);
[eL2v,eH1v]=errorFEM3D({v,vx,vy,vz},uh(Ndof+1:2*Ndof,end),T,k);
[eL2w,eH1w]=errorFEM3D({w,wx,wy,wz},uh(2*Ndof+1:end,end),T,k);

eL2 = eL2u+eL2v+eL2w;
eH1 = eH1u+eH1v+eH1w;

if ~exist('external')
    disp(eL2);
    disp(eH1);
else
    errors = [errors; eL2 eH1];
end