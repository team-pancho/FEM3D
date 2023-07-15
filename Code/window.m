function [wf,wfp,wfpp]=window(a,b,c,d)
%
% [wf,wfp,wfpp]=window(a,b,c,d)
%
% Input:
% a,b,c,d       : Numbers such that
%                 wf(x)=0 on (-infty,a]
%                 wf(x)=inc on [a,b]
%                 wf(x)=1 on [b,c]
%                 wf(x)=dec on [c,d]
%                 wf(x)=0 on [d, infty]
%
% Output:
% wf            : A function handle
% wfp, wfpp     : Function handles (derivatives of wf)
%
% Last modified: November 21, 2017

HH =  @(x) x.^5.*(1-5*(x-1)+15*(x-1).^2-35*(x-1).^3+70*(x-1).^4-...
    126*(x-1).^5).*(x>0).*(x<1)+(x>=1);
HHp = @(x) (5*x.^4.*...
    (1-5*(x-1)+15*(x-1).^2-35*(x-1).^3+70*(x-1).^4-126*(x-1).^5)...
    +x.^5.*(-5+30*(x-1)-105*(x-1).^2+280*(x-1).^3-630*(x-1).^4))...
    .*(x>0).*(x<1);
HHpp= @(x) -1260*(x-1).^4.*(9*x-4).*(x.^3).*(x>0).*(x<1);
wf  = @(x) HH((x-a)./(b-a)).*HH((d-x)./(d-c));
wfp = @(x) 1/(b-a)*HHp((x-a)/(b-a)).*HH((d-x)./(d-c))...
    -1/(d-c)*HH((x-a)/(b-a)).*HHp((d-x)./(d-c));
wfpp = @(x) 1/(b-a)^2*HHpp((x-a)/(b-a)).*HH((d-x)./(d-c))...
    -2/(b-a)*1/(d-c)*HHp((x-a)/(b-a)).*HHp((d-x)./(d-c))...
    +1/(d-c)^2*HH((x-a)/(b-a)).*HHpp((d-x)./(d-c));