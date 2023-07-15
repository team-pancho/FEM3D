function [ S ] = meshTetRefinement( T )

% function [ S ] = meshTetRefinement( T )
% 
% Input :
%        T : a basic Tetrahedral mesh
% 
% Output :
%        S : a new mesh where each element has been refined into 8 new
%            tetrahedra
% 
% Last Modified: May 20, 2016

Nvert = size(T.coordinates,2);
Nelt = size(T.elements,2);
NDir = size(T.dirichlet,2);
NNeu = size(T.neumann,2);
T = edgesAndFaces(T);

S.coordinates = [T.coordinates   (T.coordinates(:,T.edges(1,:))...
                                  + T.coordinates(:, T.edges(2,:))) /2];
newVertoldElt = abs(T.edgebyelt) + Nvert;
allVert = [T.elements; newVertoldElt];

% following pattern taken from Liu and Joe '96, Quality Local Refinement of
% Tetrahedral Meshes Based on 8-subtetrahedron subdivision
newEltPattern = [1 5 7 10 5 2 6 8 7 6 3 9 10 8 9 4 ...% these are good
                 5 7 10 8 5 8 6 7 9 7 6 8 9 10 7 8]'; % these can be changed 
                                                      % based on
                                                      % deformation of the
                                                      % tets (maybe even
                                                      % element by element)
S.elements = reshape(allVert(newEltPattern,:),4, 8*Nelt);

faceByEdges = abs(T.edgebyelt([1 2 3 1 4 6 3 5 6 4 2 5],:));
faceByEdges = reshape(faceByEdges,3,4*Nelt) + Nvert;
x = T.facebyelt(:);
orients = T.faceorient(:);

rot = [1 1 3 3 2 2;...
       2 3 1 2 3 1;...
       3 2 2 1 1 3];
if NDir
    dir = find(T.faces(4,:) ==1);
    dirPos = find(ismember(x,dir) ==1);
    dirFaceCopy = T.faces;
    for i = 1:length(dirPos)
        dirFaceCopy(rot(:, orients(dirPos(i))), x(dirPos(i)))...
            = dirFaceCopy(1:3, x(dirPos(i)));
    end

    dirInElts = dirFaceCopy(1:3, x(dirPos));
    newDir = [dirInElts(1,:); faceByEdges(1,dirPos); ...
              dirInElts(2,:); faceByEdges(2,dirPos); ...
              dirInElts(3,:); faceByEdges(3,dirPos)];
    newDir = newDir([1 2 6 2 3 4 6 4 5 2 4 6],:);
    newDir = reshape(newDir,3,4*length(dirPos));

    for i = 1:length(dirPos)
        S.dirichlet(:,4*(i-1) + 1:4*i) = ...
            newDir(rot(:, orients(dirPos(i))), 4*(i-1) + 1:4*i);
    end
else 
    S.dirichlet = [];
end

if NNeu
    neu = find(T.faces(4,:) == 2);
    neuPos = find(ismember(x,neu) ==1);
    neuFaceCopy = T.faces;
    for i = 1:length(neuPos)
        neuFaceCopy(rot(:, orients(neuPos(i))), x(neuPos(i)))...
            = neuFaceCopy(1:3, x(neuPos(i)));
    end

    neuInElts = neuFaceCopy(1:3, x(neuPos));
    newNeu = [neuInElts(1,:); faceByEdges(1,neuPos); ...
              neuInElts(2,:); faceByEdges(2,neuPos); ...
              neuInElts(3,:); faceByEdges(3,neuPos)];
    newNeu = newNeu([1 2 6 2 3 4 6 4 5 2 4 6],:);
    newNeu = reshape(newNeu,3,4*length(neuPos));

    for i = 1:length(neuPos)
        S.neumann(:,4*(i-1) + 1:4*i) = ...
            newNeu(rot(:, orients(neuPos(i))), 4*(i-1) + 1:4*i);
    end
else
    S.neumann = [];
end

S.bdtag = T.bdtag;

% to be changed...
% AllBDFaces = [S.dirichlet, S.neumann];
% AllFaceCoords = S.coordinates(:,AllBDFaces);
% AllFaceCoords = reshape(AllFaceCoords, 3,3, size(AllBDFaces,2));
% S.bdface{1} = AllBDFaces(:,sum(AllFaceCoords(3,:,:))==0);
% S.bdface{2} = AllBDFaces(:,sum(AllFaceCoords(3,:,:))==3);
% S.bdface{3} = AllBDFaces(:,sum(AllFaceCoords(2,:,:))==0);
% S.bdface{4} = AllBDFaces(:,sum(AllFaceCoords(2,:,:))==3);
% S.bdface{5} = AllBDFaces(:,sum(AllFaceCoords(1,:,:))==0);
% S.bdface{6} = AllBDFaces(:,sum(AllFaceCoords(1,:,:))==3);

end

