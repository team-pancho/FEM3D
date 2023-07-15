function fh=loadVector3D(f,T,k)

% fh=loadVector3D(f,T,k)
%
% Input:
%    f  : vectorized function of three variables
%                       or
%         array of vectorized functions of three vars {f1,f2,f3,...}
%    T  : enhanced triangulation
%    k  : polynomial degree
% Output:
%    fh : dim P_k(T_h) load vector
%                  or matrix of load vectors
% Last modified: July 16, 2015

formula=quadratureFEM(3*k,3);
x=T.coordinates(1,:); x=formula(:,1:4)*x(T.elements);
y=T.coordinates(2,:); y=formula(:,1:4)*y(T.elements);
z=T.coordinates(3,:); z=formula(:,1:4)*z(T.elements);
Pw=bsxfun(@times,formula(:,5),bernstein3D(formula(:,2),formula(:,3),formula(:,4),k));
loc2glob=DOF3D(T,k);
dimVh=max(loc2glob(:));
if iscell(f)
    fh=zeros(dimVh,length(f));
    for i=1:length(f)
        fhi=bsxfun(@times,T.volume,Pw'*f{i}(x,y,z));
        fh(:,i)=accumarray(loc2glob(:),fhi(:));
    end
else
    fh=bsxfun(@times,T.volume,Pw'*f(x,y,z));% dk x Nelt matrix with local tests
    fh=accumarray(loc2glob(:),fh(:));       
end
return


