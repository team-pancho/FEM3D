function [localDOF,Cols]=DOF3D(T,k)

% function [localDOF,Cols]=DOF3D(T,k)
% Input:
%      T        : Data structure, enhanced FEM tetrahedrization
% Output:
%      localDOF : dimPk x Nelt matrix with local DOF
%      Cols     : 3D array for assembly
% Last modified: May 1, 2015  

Nnod=size(T.coordinates,2);
Nelt=size(T.elements,2);
Nedg=size(T.edges,2);
Nfac=size(T.faces,2);

fulldim=((k+3)*(k+2)*(k+1))/6;
dimEdges=(k-1)*(k>1);
dimFaces=((k-1)*(k-2))/2*(k>2);
dimInt  =((k-1)*(k-2)*(k-3))/6*(k>3);

vDOF=T.elements;
eDOF=[];
if k>1
    eDOF=(abs(T.edgebyelt)-1)*dimEdges+Nnod;
    eDOF=eDOF(:)';
    eDOF=bsxfun(@plus,(1:dimEdges)',eDOF);
    negative=find(sign(T.edgebyelt(:))==-1);
    eDOF(:,negative)=eDOF(dimEdges:-1:1,negative);
    eDOF=reshape(eDOF,[6*dimEdges,Nelt]);
end
fDOF=[];
if k>2
    fDOF=(T.facebyelt-1)*dimFaces+Nnod+dimEdges*Nedg;
    fDOF=fDOF(:)';
    fDOF=bsxfun(@plus,(1:dimFaces)',fDOF);
    rotation=rotateFace(k);
    for j=2:6
        perm=find(T.faceorient==j);
        fDOF(rotation(:,j),perm)=fDOF(:,perm);
    end
    fDOF=reshape(fDOF,[4*dimFaces,Nelt]);
end
tDOF=[];
if k>3
    tDOF=((1:Nelt)-1)*dimInt+Nnod+dimEdges*Nedg+dimFaces*Nfac;
    tDOF=bsxfun(@plus,(1:dimInt)',tDOF);
end
localDOF=[vDOF;eDOF;fDOF;tDOF];
Cols=repmat(localDOF(:)',[fulldim 1]);
Cols=reshape(Cols,[fulldim,fulldim,Nelt]);

return

function rotation=rotateFace(k)

index = 1;
for ind2 = 1:k-2
    for ind3 = 1:k-1-ind2
        ind1 = k-ind2-ind3;
        table(index,:) = [ind1 ind2 ind3];
        index = index +1;
    end
end
perm = [ 1 2 3;...
         1 3 2;...
         3 1 2;...
         3 2 1;...
         2 3 1;...
         2 1 3];
dim=((k-1)*(k-2))/2;
rotation=zeros(dim,6);
for j= 1:6
    [~,list] = ismember(table,table(:,perm(j,:)),'rows');
    rotation(:,j) = list;
end
return

