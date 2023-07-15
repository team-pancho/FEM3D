function [fh,M,dof,ML]=bdWork3D(f,T,k,list,goal,testmode,lambda)

% [fh,M,dof]=bdWork3D(f,T,k,list,goal,tstmod,lambda)
%
% Input: 
%     f      : cell array with vectorized functions
%     T      : enhanced tetrahedrization
%     k      : polynomial degree
%     list   : row vector (sublist of the set of faces)
%     goal   : 'test', 'proj', 'eval'
%     tstmod : 'sc' (scalar), 'fl' (flux), 'tr' (traction)
%     lambda : array of vectorized functions of three variables 
%              of one vectorized function of three variables
%              ALTERNATIVELY, any non-functional input (treated as empty)
% Output:
%     fh     : dim V_h long vector (see later to see what it does)
%     M      : dim V_h x dim V_h matrix 
%                 int_{\Gamma_list} P_i P_j,  
%                 \Gamma_list=\cup_{i\in list} F_i 
%     dof    : dim P_k(F) x #list matrix with DOF for faces F_i, i\in list
%     ML      : dim V_h x dim V_h matrix 
%                 int_{\Gamma_list} lambda P_i P_j, 
%                 of int_{\Gamma_list} (lambda\dot n) P_i P_j
%
% Last modified: January 19, 2017

% Quadrature points and geometric quantities

form = quadratureFEM(3*k,2);

x = T.coordinates(1,:);   
y = T.coordinates(2,:);   
z = T.coordinates(3,:); 
x = form(:,1:3)*x(T.faces(1:3,list));  
y = form(:,1:3)*y(T.faces(1:3,list));
z = form(:,1:3)*z(T.faces(1:3,list));

P = bernstein2D(form(:,2),form(:,3),k);
Pw = bsxfun(@times,form(:,4),P); 

normx = T.normals(1,list); 
normy = T.normals(2,list);
normz = T.normals(3,list);
areas = T.area(list);

[dof,Assem]=computeBDDOF3D(T,k,list);
dimVh=dimFEMspace(T,k);

% Testing on separate faces

if testmode=='sc'
    for ell=1:length(f)
        fh{ell}=bsxfun(@times,areas,f{ell}(x,y,z));
    end
elseif testmode=='fl'
    fh{1}=bsxfun(@times,normx,f{1}(x,y,z))...
            +bsxfun(@times,normy,f{2}(x,y,z))...
            +bsxfun(@times,normz,f{3}(x,y,z));
elseif testmode=='tr'
    fh{1} = bsxfun(@times,normx,f{1}(x,y,z))...
                +bsxfun(@times,normy,f{2}(x,y,z))...
                +bsxfun(@times,normz,f{3}(x,y,z));
    fh{2} = bsxfun(@times,normx,f{2}(x,y,z))...
                +bsxfun(@times,normy,f{4}(x,y,z))...
                +bsxfun(@times,normz,f{5}(x,y,z));
    fh{3} = bsxfun(@times,normx,f{3}(x,y,z))...
                +bsxfun(@times,normy,f{5}(x,y,z))...
                +bsxfun(@times,normz,f{6}(x,y,z));
end
for ell=1:length(fh)
    fh{ell}=Pw'*fh{ell};
    fh{ell}=accumarray(dof(:),fh{ell}(:),[dimVh,1]);
end
if testmode=='sc'
    FF=sparse(dimVh,length(f));
    for ell=1:length(f);
        FF(:,ell)=fh{ell};
    end
    fh=FF;
elseif testmode=='fl'
    fh=fh{1};
elseif testmode=='tr'
    fh=[fh{1} fh{2} fh{3}];
end

% Boundary mass matrix with unit coefficient

M=kron(areas,P'*Pw);
Rows=permute(Assem,[2 1 3]);
M=sparse(Rows(:),Assem(:),M(:),dimVh,dimVh);


% L2 projection: only used when goal='proj' mode

if goal=='proj'
    active=unique(dof(:));
    fh(active,:)=M(active,active)\fh(active,:);
end

% Boundary mass matrix with variable coefficient

if iscell(lambda) || isa(lambda,'function_handle')
    if iscell(lambda)
        Lambda=bsxfun(@times,normx,lambda{1}(x,y,z))...
               +bsxfun(@times,normy,lambda{2}(x,y,z))...
               +bsxfun(@times,normz,lambda{3}(x,y,z));
    else
        Lambda=bsxfun(@times,areas,lambda(x,y,z));
    end
    dk=size(P,2);
    ML=zeros(dk,dk*length(list));
    for q=1:size(form,1)
        ML=ML+kron(Lambda(q,:),form(q,4)*P(q,:)'*P(q,:));
    end
    ML=sparse(Rows(:),Assem(:),ML(:),dimVh,dimVh);
else
    ML = 0;
end  

return

