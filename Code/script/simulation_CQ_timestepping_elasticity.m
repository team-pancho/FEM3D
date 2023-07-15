% A simulation to couple CQ time stepping with 3D FEM for elasticity and
% make a movie
% Last modified: May 20, 2016

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

%plane wave parameters
d=[1 1 1]/sqrt(3);
e=d;
a=1;
tlag=0.3;

% material parameters & exact solution
mu  =@(x,y,z) 1+0*x;
lam =@(x,y,z) 1+0*x;
rho =@(x,y,z) 1+0*x;

lamx =@(x,y,z) 0*x;
lamy =@(x,y,z) 0*y;
lamz =@(x,y,z) 0*z;

mux =@(x,y,z) 0*x;
muy =@(x,y,z) 0*y;
muz =@(x,y,z) 0*z;

% wavespeed
if e==d
    c=sqrt((lam(1,1,1)+2*mu(1,1,1))/rho(1,1,1));
else
    c=sqrt(mu(1,1,1)/rho(1,1,1));
end


% u
u=@(x,y,z,t)  a*e(1)*signal(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
ux=@(x,y,z,t) -d(1)*a*e(1)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uy=@(x,y,z,t) -d(2)*a*e(1)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uz=@(x,y,z,t) -d(3)*a*e(1)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uxx=@(x,y,z,t) d(1)^2*a*e(1)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uyy=@(x,y,z,t) d(2)^2*a*e(1)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uzz=@(x,y,z,t) d(3)^2*a*e(1)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uxy=@(x,y,z,t) d(1)*d(2)*a*e(1)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uxz=@(x,y,z,t) d(1)*d(3)*a*e(1)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uyz=@(x,y,z,t) d(2)*d(3)*a*e(1)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
ut =@(x,y,z,t) a*e(1)*c*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
utt=@(x,y,z,t) a*e(1)*c^2*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));

% v
v=@(x,y,z,t)  a*e(2)*signal(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vx=@(x,y,z,t) -d(1)*a*e(2)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vy=@(x,y,z,t) -d(2)*a*e(2)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vz=@(x,y,z,t) -d(3)*a*e(2)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vxx=@(x,y,z,t) d(1)^2*a*e(2)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vyy=@(x,y,z,t) d(2)^2*a*e(2)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vzz=@(x,y,z,t) d(3)^2*a*e(2)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vxy=@(x,y,z,t) d(1)*d(2)*a*e(2)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vxz=@(x,y,z,t) d(1)*d(3)*a*e(2)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vyz=@(x,y,z,t) d(2)*d(3)*a*e(2)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vt =@(x,y,z,t) a*e(2)*c*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vtt=@(x,y,z,t) a*e(2)*c^2*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));

% w
w=@(x,y,z,t)  a*e(3)*signal(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wx=@(x,y,z,t) -d(1)*a*e(3)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wy=@(x,y,z,t) -d(2)*a*e(3)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wz=@(x,y,z,t) -d(3)*a*e(3)*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wxx=@(x,y,z,t) d(1)^2*a*e(3)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wyy=@(x,y,z,t) d(2)^2*a*e(3)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wzz=@(x,y,z,t) d(3)^2*a*e(3)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wxy=@(x,y,z,t) d(1)*d(2)*a*e(3)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wxz=@(x,y,z,t) d(1)*d(3)*a*e(3)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wyz=@(x,y,z,t) d(2)*d(3)*a*e(3)*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wt =@(x,y,z,t) a*e(3)*c*signalp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wtt=@(x,y,z,t) a*e(3)*c^2*signalpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));

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

uh=real(ifft(uh,[],2)); uh=bsxfun(@times,uh,R.^(-0:M));        % v

%file stuff
DateAndTime=datetime('now');
TimeStamp=datestr(datenum(DateAndTime),'yy-mm-dd-HH-MM');
FileStr=...
    strcat('/Users/allanhungria/Dropbox/HDG3dHelmholtz/scripts/movies/',...
            TimeStamp,'_ElasticHelmholtzMovie');
ntime_=0:size(uh,2)-1;

close all

% PlotTransientElasticity(uh,T,k,varargin)


if k==1
    time_=kappa*ntime_;
    Muh=max(max(uh));
    muh=min(min(uh));
    Mte=max(max(T.coordinates));
    mte=min(min(T.coordinates));
    LIM=[min([muh mte muh+mte]),max([Muh Mte Muh+Mte])];
    tn=0;
    for t = time_
        tn=tn+1;
        nowDOF=uh(:,tn);
        tetramesh(T.elements',T.coordinates'...
                  +reshape(nowDOF,length(nowDOF)/3,3))
        xlim(LIM)
        ylim(LIM)
        zlim(LIM)
        grid on
        pause(0.1)
    end
end