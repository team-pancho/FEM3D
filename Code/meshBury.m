function T = meshBury(T,L)

% T=meshBury(T,L)
%
% Input:
%    T   : tetrahedral mesh
%    L   : level
% Output:
%    T   : tetrahedral mesh - boundary face fields have been subdivided 
%          based on the location of the z-coordinates of the barycenters 
%          with respect to L
% Last modified: July 8, 2016

Facenum = size(T.bdface,2);

for i=1:Facenum 
    FaceCoords = T.coordinates(:,T.bdface{i});
    n = size(FaceCoords,2);

    F1 = FaceCoords(:,1:3:n);
    F2 = FaceCoords(:,2:3:n);
    F3 = FaceCoords(:,3:3:n);

    BaryC = (F1+F2+F3)/3;

    M = find(BaryC(3,:)<=L);

    T.bdface{i+Facenum} = T.bdface{i}(:,M);
    T.bdface{i}(:,M)=[];
    T.bdtag{i+Facenum}=[ 'below ' T.bdtag{i}];
    T.bdtag{i}=['above ' T.bdtag{i}];
end

T.dirichlet=[];
T.neumann=[];
for i=Facenum+1:2*Facenum
    T.dirichlet=[T.dirichlet T.bdface{i}];
end
for i=1:Facenum
    T.neumann=[T.neumann T.bdface{i}];
end

empty=[];
for i=1:2*Facenum
    if size(T.bdface{i},1) == 0 || size(T.bdface{i},2) == 0 
        if  (i+Facenum)<=size(T.bdface,2) && size(T.bdface{i+Facenum},1)>=1
            T.bdface([i,i+Facenum]) = T.bdface([i+Facenum,i]);
            T.bdtag([i,i+Facenum]) = T.bdtag([i+Facenum,i]);
            empty = [empty i+Facenum];
        else
            empty = [empty i];
        end
    end
end
T.bdface(empty)=[];
T.bdtag(empty)=[];

return


