function T=edgesAndFaces(T)

% T=edgesAndFaces(T)
% Input: 
%     T    : basic tetrahedrization data structure (3D tets)
% Output:
%     T    : geometric fields added 
%            T.edges, T.edgebyelt, T.faces,T.faceorient,T.facebyelt
%
% Last modified: May 4, 2015.

N=size(T.coordinates,2);
Nelt=size(T.elements,2);
Ndir = size(T.dirichlet,2);
Nneu = size(T.neumann,2);


if Ndir
    edgesDir = sparse(T.dirichlet,T.dirichlet([2 3 1],:),1); % all dir face edges
    [idir, jdir] = find(triu(edgesDir+edgesDir'));
    diredges=[idir';jdir'];
    edgesDir=sparse(diredges(1,:),diredges(2,:),1,N,N);
else
    edgesDir=sparse(1,1,0,N,N);
    idir=[];
    jdir=[];   
end


if Nneu
    edgesNeu = sparse(T.neumann,T.neumann([2 3 1],:),1); % all neu face edges
    [ineu, jneu] = find(triu(edgesNeu+edgesNeu')==2);
    neuedges=[ineu';jneu'];
    edgesNeu=sparse(neuedges(1,:),neuedges(2,:),1,N,N);
else 
    edgesNeu=sparse(1,1,0,N,N);
    ineu = [];
    jneu = [];
end

    
edgesAll=sparse(T.elements([1 2 1 2 3 1],:),T.elements([2 3 3 4 4 4],:),1,N,N);
edgesAll=sign(triu(edgesAll+edgesAll'));
edgesInt=edgesAll-edgesDir-edgesNeu;

[iint, jint]=find(edgesInt);
T.edges=[iint jint zeros(size(iint,1),1);...
         idir jdir ones(size(idir,1),1);...
         ineu jneu 2*ones(size(ineu,1),1)]';
Nedge=size(T.edges,2);
edgesAll=sparse(T.edges(1,:),T.edges(2,:),1:Nedge,N,N);
edgesAll=edgesAll-edgesAll';

T.edgebyelt=...
    [full(edgesAll(sub2ind([N,N],T.elements(1,:),T.elements(2,:))));...
     full(edgesAll(sub2ind([N,N],T.elements(2,:),T.elements(3,:))));...
     full(edgesAll(sub2ind([N,N],T.elements(1,:),T.elements(3,:))));...
     full(edgesAll(sub2ind([N,N],T.elements(2,:),T.elements(4,:))));...
     full(edgesAll(sub2ind([N,N],T.elements(3,:),T.elements(4,:))));...
     full(edgesAll(sub2ind([N,N],T.elements(1,:),T.elements(4,:))))];

 
shape=[1 1 1 4;...
       2 2 3 2;...
       3 4 4 3];
faces=reshape(T.elements(shape(:),:),3,4*Nelt);
faces=sort(faces,1)';
[allfaces,~,j]=unique(faces,'rows');

% Lists of interior and boundary faces with references
bdfaces=sort([T.dirichlet T.neumann],1)';
[intfaces,i]=setdiff(allfaces,bdfaces,'rows');
[bdfaces,ii,jj]=intersect(allfaces,bdfaces,'rows');

%sizes
nintfaces=size(intfaces,1);
nbdfaces =size(bdfaces,1);
nfaces   =nintfaces+nbdfaces;

T.faces=[intfaces' T.dirichlet T.neumann;
         zeros(1,nintfaces) ones(1,Ndir) 2*ones(1,Nneu)];

% Backward referencing to construct T.facebyelt
u=zeros(nfaces,1);
v=nintfaces+1:nintfaces+nbdfaces;
u(i)=1:nintfaces;
u(ii)=v(jj);
j=u(j);         % reconstructing list of all faces (4*Nelt) now with global numbering 
T.facebyelt=reshape(j,[4 Nelt]);


% Matrix with permutation order
eq = @(u,v) sum(u==v,1)==3; % checks what rows are equal
rot=[1 1 3 3 2 2;...
     2 3 1 2 3 1;...
     3 2 2 1 1 3];       %permutations

T.faceorient=zeros(4,Nelt);
for f=1:4    % counter over faces
    faceGlobal=T.faces(1:3,T.facebyelt(f,:));
    faceLocal =T.elements(shape(:,f),:);
    for j=1:6
        T.faceorient(f,:)=T.faceorient(f,:)+j*eq(faceGlobal,faceLocal(rot(:,j),:));
    end
end

return