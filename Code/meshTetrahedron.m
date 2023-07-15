function T=meshTetrahedron(N,varargin)

% T=reftetMesh(N)
% T=reftetMesh(N,listDir)
% input
%    N  : number of partitions in the x, y and z directions
%    listDir : list of Dirichlet faces
% output
%    T   : tetrahedral partition data structure
%          default: all Dirichlet
%
% last modified: March 4, 2016

TQ=meshCube(N,N,N);
ep=0.5/N;
stay=find(sum(TQ.coordinates)<=1+ep);
nvert=size(TQ.coordinates,2);
nvertTet=length(stay);
corresp=zeros(nvert,1);
corresp(stay)=1:nvertTet;
todelete=@(A) find(prod(A)==0);

% Coordinates and elements

T.coordinates=TQ.coordinates(:,stay);
T.elements   =corresp(TQ.elements);
T.elements(:,todelete(T.elements))=[];

% Coordinate faces

T.bdface{1}=corresp(TQ.bdface{3});
T.bdface{2}=corresp(TQ.bdface{1});
T.bdface{3}=corresp(TQ.bdface{5});
T.bdface{1}(:,todelete(T.bdface{1}))=[];
T.bdface{2}(:,todelete(T.bdface{2}))=[];
T.bdface{3}(:,todelete(T.bdface{3}))=[];

T.bdtag{1}='left (y=0)';
T.bdtag{2}='bottom (z=0)';
T.bdtag{3}='back (x=0)';

% Lid

ep=0.05/N;
notonlid=find(abs(sum(T.coordinates)-1)>ep);
corresp=(1:nvertTet)';
corresp(notonlid)=0;
group1=corresp(T.elements([1 2 3],:));
group2=corresp(T.elements([2 3 4],:));
group3=corresp(T.elements([3 4 1],:));
group4=corresp(T.elements([4 1 2],:));
group1(:,todelete(group1))=[];
group2(:,todelete(group2))=[];
group3(:,todelete(group3))=[];
group4(:,todelete(group4))=[];
T.bdface{4}=[group1 group2 group3 group4];
T.bdtag{4}='lid (x+y+z=1)';

% Dirichlet/Neumann faces

if nargin==1
    listDir=[1 2 3 4];
    listNeu=[];
else
    listDir=varargin{1};
    listNeu=1:4;
    listNeu(listDir)=[];
end
T.dirichlet=[];
T.neumann=[];
for i=listDir
    T.dirichlet=[T.dirichlet T.bdface{i}];
end
for i=listNeu
    T.neumann=[T.neumann T.bdface{i}];
end

end

