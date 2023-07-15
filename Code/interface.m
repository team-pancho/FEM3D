function T = interface(T)

%function T = interface(T)
% 
% Input:
%       T: data structure of a subtetrahedrization of a larger 
%          tetrahedrization 
% 
% Output:
%       T: updated data structure with the following changes:
%          (a) a new field T.interface which identifies false "boundary" faces
%          (b) an option of three in the last row in T.faces of the above
%              faces 
%          (c) T.normals has been updated so that all normal vectors on the
%              interface point in the same direction
% 
% Last Modified: October 21, 2016

freq = accumarray(T.facebyelt(:),1)';
nonintface = find(freq==1);
nonbdface = find(T.faces(4,:)==0);
interfaces = intersect(nonintface,nonbdface);
T.interface = T.faces(1:3,interfaces);
T.faces(4,interfaces) = 3;

A=T.faceorient;
for row=[1 3]
    for perm=[2 4 6]
        A(row,A(row,:)==perm)=10; % 10 is fine, 0 is not fine
    end
    A(row,A(row,:)~=10)=0;
end
for row=[2 4]
    for perm=[1 3 5]
        A(row,A(row,:)==perm)=10; % 10 is fine, 0 is not fine
    end
    A(row,A(row,:)~=10)=0;
end
A(A==10)=1;

[~,WhichFaceIsIn]=ismember(T.facebyelt,interfaces);
Wrong=WhichFaceIsIn.*(1-A); % marks is the face is in and has the wrong orientation
Wrong=unique(Wrong(:)); 
if Wrong(1)==0; 
    Wrong=Wrong(2:end); 
end
change=interfaces(Wrong);
T.normals(:,change) = -T.normals(:,change);
end
