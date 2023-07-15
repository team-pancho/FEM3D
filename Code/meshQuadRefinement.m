function [Qnew]=meshQuadRefinement(Q)

% [Qnew]=meshQuadRefinement(Q)
%
% Input: 
%   Q      : basic quad partition
% Output:
%   Qnew     : uniformly refined quad partition
% 
% Last Modified: July 1, 2016

Ncubes=size(Q.elements,2);
Nnodes=size(Q.coordinates,2);

% Edge connectivity

edges=[1 2;2 3;3 4;4 1;...
       5 6;6 7;7 8;8 5;... 
       1 5; 2 6;3 7; 4 8]';
econ=reshape(Q.elements(edges(:),:),2,12*Ncubes);
econ=sort(econ);
whichedge=zeros(12,Ncubes);
[econ,i,whichedge(:)]=unique(econ','rows');  % econ = Nedges x 2
% whichedge=whichedge';     % whichedge(e,:) = edges for element #e
Nedges=size(econ,1);
econ = econ';

% Face connectivity 

faces=[1 2 3 4;...
       5 8 7 6;...
       1 5 6 2;...
       4 3 7 8;...
       1 4 8 5;...
       2 6 7 3]';
fcon=reshape(Q.elements(faces(:),:),4,6*Ncubes);
ffcon=fcon';
fcon=sort(fcon);
whichface=zeros(6,Ncubes);
[fcon,i,whichface(:)]=unique(fcon','rows');  % fcon = Nfaces x 4
% whichface=whichface';     % whichface(e,:) = faces for element #e
Nfaces=size(fcon,1);   
ffcon=ffcon(i,[2 4]);
fcon = fcon';


% Coordinates of new points

Qnew.coordinates=[Q.coordinates ...
   ( Q.coordinates(:,econ(1,:))+Q.coordinates(:,econ(2,:)))/2 ...           % edges
   ( Q.coordinates(:,fcon(1,:))+Q.coordinates(:,fcon(2,:))...               % faces
    +Q.coordinates(:,fcon(3,:))+Q.coordinates(:,fcon(4,:)))/4 ...
   ( Q.coordinates(:,Q.elements(1,:))+Q.coordinates(:,Q.elements(2,:))...   % elements
    +Q.coordinates(:,Q.elements(3,:))+Q.coordinates(:,Q.elements(4,:))...
    +Q.coordinates(:,Q.elements(5,:))+Q.coordinates(:,Q.elements(6,:))...
    +Q.coordinates(:,Q.elements(7,:))+Q.coordinates(:,Q.elements(8,:)))/8];

% Element subdivision

elt=[Q.elements;...
        Nnodes+whichedge;...
        Nnodes+Nedges+whichface;...
        Nnodes+Nedges+Nfaces+(1:Ncubes)];
pattern=[1 9 21 12 17 23 27 25;...
    9 2 10 21 23 18 26 27;...
    12 21 11 4 25 27 24 20;...
    21 10 3 11 27 26 19 24;...
    17 23 27 25 5 13 22 16;...
    25 27 24 20 16 22 15 8;...
    23 18 26 27 13 6 14 22;...
    27 26 19 24 22 14 7 15]';
Qnew.elements=reshape(elt(pattern,:),8,8*Ncubes);

% PARTITION OF BOUNDARY FACES

edges=[1 2 3 4;...
       2 3 4 1];
pattern=[1 5 9 8;...
    5 2 6 9;...
    8 9 7 4;...
    9 6 3 7];

% Location of edges and faces of Dirichlet elements

Ndir=size(Q.dirichlet,2); 
if Ndir>0
    listDir=reshape(Q.dirichlet(edges(:),:),2,4*Ndir);
    listDir=sort(listDir',2);       % edges of Dirichlet faces (with repetitions)
    [listDir,i,j]=unique(listDir,'rows');
    [aux,ii,jj]=intersect(econ',listDir,'rows');
    whichDiredge=reshape(ii(j),4,Ndir);
    [aux,i,j]=intersect(fcon',sort(Q.dirichlet',2),'rows');
    whichDirface(j) = i;
    elt=[Q.dirichlet;...
        Nnodes+whichDiredge;...
        Nnodes+Nedges+whichDirface];
    Qnew.dirichlet=reshape(elt(pattern',:),4,4*Ndir);
else
    Qnew.dirichlet=[];
end

% Location of edges and faces of Neumann elements

Nneu=size(Q.neumann,2); 
if Nneu>0
    listNeu=reshape(Q.neumann(edges(:),:),2,4*Nneu);
    listNeu=sort(listNeu',2);       % edges of Neumann faces (with repetitions)
    [listNeu,i,j]=unique(listNeu,'rows');
    [aux,ii,jj]=intersect(econ',listNeu,'rows');
    whichNeuedge=reshape(ii(j),4,Nneu);
    [aux,i,j]=intersect(fcon',sort(Q.neumann',2),'rows');
    whichNeuface(j)=i;
    elt=[Q.neumann;...
        Nnodes+whichNeuedge;...
        Nnodes+Nedges+whichNeuface];
    Qnew.neumann=reshape(elt(pattern',:),4,4*Nneu);
else
    Qnew.neumann=[];
end

if isfield(Q, 'bdface')
    disp('WARNING: at this time partition.m is not equipped to update');
    disp('the fields bdface and bdtag.  As such, they no longer exist.');
end    
return


