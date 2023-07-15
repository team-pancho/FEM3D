function u = CQoperator(F,g,k,dim,varargin)
% u = CQoperator(F,g,k,dim)
% u = CQoperator(F,g,k,dim,p)
%
% Input:
%    F : dim x 1 valued transfer function. Function handle of
%        s : complex number, X : d2 x 1 complex valued vector
%    g : Time values of input (d2 x (M+1) matrix)
%    k : Time-step
%  dim : Dimension of output for F
%    p : Transfer function of multistep method. Function handle of s
%
% Output:
%    u : dim x (M+1) matrix corresponding to F(\partial_k) g           
%
% Last Modified : September 19, 2018

if nargin==4
    p = @(z) 2*(1-z)./(1+z);     % The default is TR 
else
    p=varargin{1};
end

M  = size(g,2)-1;              % M = number of time-steps
omega = exp(2*pi*1i/(M+1));
R = eps^(0.5/(M+1));

h = bsxfun(@times,g,R.^(0:M));    
h = fft(h,[],2);               % DFT by columns (\hat H) 
u = zeros(dim,M+1);
parfor l=0:floor((M+1)/2)
    u(:,l+1)=F(p(R*omega^(-l))/k,h(:,l+1));%#ok     % \hat v   
end
u(:,M+2-(1:floor(M/2)))=conj(u(:,2:floor(M/2)+1));
u=real(ifft(u,[],2));                           % v
u=bsxfun(@times,u,R.^(-(0:M)));

return
