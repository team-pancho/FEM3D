function [u,v,w,ut,vt,wt,utt,vtt,wtt,...
                ux,uy,uz,vx,vy,vz,wx,wy,wz,...
                uxx,uxy,uxz,uyy,uyz,uzz,...
                vxx,vxy,vxz,vyy,vyz,vzz,...
                wxx,wxy,wxz,wyy,wyz,wzz] ...
                = vectorPlaneWave(a,d,c,tlag,signal)

%        [u,v,w,ut,vt,wt,utt,vtt,wtt,...
%               ux,uy,uz,vx,vy,vz,wx,wy,wz,...
%               uxx,uxy,uxz,uyy,uyz,uzz,...
%               vxx,vxy,vxz,vyy,vyz,vzz,...
%               wxx,wxy,wxz,wyy,wyz,wzz]    ...
%               = vectorPlaneWave(a,d,c,tlag,signal)
%
% Input:
%     a       : vector amplitude of wave (3-vector)
%     d       : unit propagation direction vector (3-vector)
%     c       : wave speed
%     tlag    : lag time
%     signal  : 3-cell array with f, f', f''
%
%  Output:
%     u,v,w                   : function handles for a plane wave
%                               (u,v,w) = f(d.x-c(t-tlag))a
%     ux,uy,...,wz            : function handles with first derivatives
%     ut,vt,wt,utt,vtt,wtt    : function handles with time derivs
%     uxx,uxy,...,wzz         : function handles with second derivatives
%                               
% Last modified May 20, 2016.

f = signal{1};
fp = signal{2};
fpp = signal{3};

% u
u=@(x,y,z,t)  a(1)*f(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
ux=@(x,y,z,t) -d(1)*a(1)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uy=@(x,y,z,t) -d(2)*a(1)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uz=@(x,y,z,t) -d(3)*a(1)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uxx=@(x,y,z,t) d(1)^2*a(1)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uyy=@(x,y,z,t) d(2)^2*a(1)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uzz=@(x,y,z,t) d(3)^2*a(1)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uxy=@(x,y,z,t) d(1)*d(2)*a(1)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uxz=@(x,y,z,t) d(1)*d(3)*a(1)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
uyz=@(x,y,z,t) d(2)*d(3)*a(1)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
ut =@(x,y,z,t) a(1)*c*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
utt=@(x,y,z,t) a(1)*c^2*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));

% v
v=@(x,y,z,t)  a(2)*f(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vx=@(x,y,z,t) -d(1)*a(2)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vy=@(x,y,z,t) -d(2)*a(2)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vz=@(x,y,z,t) -d(3)*a(2)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vxx=@(x,y,z,t) d(1)^2*a(2)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vyy=@(x,y,z,t) d(2)^2*a(2)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vzz=@(x,y,z,t) d(3)^2*a(2)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vxy=@(x,y,z,t) d(1)*d(2)*a(2)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vxz=@(x,y,z,t) d(1)*d(3)*a(2)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vyz=@(x,y,z,t) d(2)*d(3)*a(2)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vt =@(x,y,z,t) a(2)*c*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
vtt=@(x,y,z,t) a(2)*c^2*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));

% w
w=@(x,y,z,t)  a(3)*f(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wx=@(x,y,z,t) -d(1)*a(3)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wy=@(x,y,z,t) -d(2)*a(3)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wz=@(x,y,z,t) -d(3)*a(3)*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wxx=@(x,y,z,t) d(1)^2*a(3)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wyy=@(x,y,z,t) d(2)^2*a(3)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wzz=@(x,y,z,t) d(3)^2*a(3)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wxy=@(x,y,z,t) d(1)*d(2)*a(3)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wxz=@(x,y,z,t) d(1)*d(3)*a(3)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wyz=@(x,y,z,t) d(2)*d(3)*a(3)*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wt =@(x,y,z,t) a(3)*c*fp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));
wtt=@(x,y,z,t) a(3)*c^2*fpp(-(x*d(1)+y*d(2)+z*d(3))+c*(t-tlag));

end

