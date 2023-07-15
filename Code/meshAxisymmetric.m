function T = meshAxisymmetric(Nx,Nz,H,r,varargin)

% T=Cylinder(Nx,Nz,H,r)
% T=Cylinder(Nx,Nz,H,r,Dir)
% Input:
%      Nx   :  8 Nx is the number of radial divisions
%      Nz   :  divisions in the z direction
%      r    :  function handle giving shape to the lateral face
%      H    :  height of cylinder in the z-direction
%      Dir  :  row vector with Dir BC (1:bottom, 2:top, 3:sides)
%              default Dir=[1 2];
% Output:
%      T    :  tetrahedrization of a cylinder
% Last modified: March 11, 2016

% Original prism

T=meshPrism(Nx,Nx,Nz,1:5);
T.coordinates(3,:)=H.*T.coordinates(3,:);
Nvert=size(T.coordinates,2);

% Symmetrized prism 

TT=T;
TT.coordinates(1:2,:)=1-TT.coordinates([2 1],:);
TT.elements=TT.elements([1 2 4 3],:);
for i=1:5
    TT.bdface{i}=TT.bdface{i}([1 3 2],:);
end

% Counter for new coordinate numbers

onInterface=unique(TT.bdface{4}(:));
notOnInterface=1:Nvert; notOnInterface(onInterface)=[];
new=1:Nvert;
new(notOnInterface)=Nvert+(1:size(notOnInterface,2));

% Merger

T.coordinates(:,new)=TT.coordinates;
T.elements=[T.elements new(TT.elements)];
T.bdface{1}=[T.bdface{1} new(TT.bdface{1})]; 
T.bdtag{1} ='bottom';
T.bdface{2}=[T.bdface{2} new(TT.bdface{2})];
T.bdtag{2} ='top';
T.bdface{3}=[T.bdface{3} T.bdface{5} new(TT.bdface{3}) new(TT.bdface{5})];
T.bdtag{3} ='all around';
T.bdface([4 5])=[];
T.bdtag([4 5])=[];

% Cylinder

T.coordinates(1:2,:)=2*T.coordinates(1:2,:)-1;   
transf = @(X,Y,Z) r(Z).*max(abs(X),abs(Y))./sqrt(X.^2+Y.^2);
rad    = transf(T.coordinates(1,:),T.coordinates(2,:),T.coordinates(3,:));
rad(isnan(rad))=0;
T.coordinates(1,:)=T.coordinates(1,:).*rad+1;
T.coordinates(2,:)=T.coordinates(2,:).*rad+1;

% Dirichlet and Neumann boundaries

if nargin==4
    listDir=[1 2];
    listNeu=3;
else
    listDir=varargin{1};
    listNeu=1:3;
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

return


