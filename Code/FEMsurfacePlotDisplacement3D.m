function FEMsurfacePlotDisplacement3D(uhtot,T,k,varargin)
% FEMsurfacePlotDisplacement3D(uhtot,T,k)
% FEMsurfacePlotDisplacement3D(uhtot,T,k,{})
%     OR     FEMsurfacePlotDisplacement3D(uhtot,T,k,movie)
% FEMsurfacePlotDisplacement3D(uhtot,T,k,movie,snapshots)
%
% Input:
%   uhtot   :      3*dimVh x Ntimes complex vector (containing uh,vh,wh)
%   T       :      Data structure. Enhanced FEM triangulation
%   k       :      Scalar. Polynomial degree
%                  is desired pass in a scalar i.e. 1.
%   movie   :      Cell element containing the following 
%                  (for default you can pass in {})
%   movie{1}:      Rotation. This can be 1 or 0, enables or disables
%                     camera rotation.
%                     Default: no rotation
%   movie{2}:      Camera angle. 1x2 vector corresponding to 
%                     [azimuth elevation] in radians.
%                     Input option: 'default'
%                     Default: [-37.5,30]
%   movie{3}:      Color for the faces in MATLAB RGB vector,
%                     Input options: 'default' OR 'flat', 
%                      OR 'grayscale', 'hotscale, ''winterscale',...
%                      OR 'green','red',...
%                      OR [1,0.3,0.8]
%                     'flat': Uniform colors taken from stress data
%                     'grayscale': Computes the grayscale of the
%                                  given coloring in the stress data
%                     'green': Uniform red color ignoring stress
%                     Default: 'flat'
%   movie{4}:      Axis on or off. This can be 1 or 0. When
%                     it is on it shows the labels too.
%                     Default: on
%   movie{5}:      Edge line style. 
%                     Input options: '-'(solid line), 'none', '--'(dashed),
%                     ':' (dotted)
%                     Default: '-'
%   movie{6}:      Number of frames per second. A positive integer.
%                     Default: 30
%   snapshots:      Enables to save snaphots.
%                   Input options: 'all',index vector
%                     'all': Saves all time steps of the snaphots of the given
%                     data
%                     index vector: Saves the snapshots at the time steps
%                     corresponding to the index vector
%                   When you save snapshots, they are all saved in the
%                   designated folder with names: snapshot_1,snapshot_2,...
%                   in png format.
% Action: plot of ...
%   displaced mesh with constant coloring
%   displaced mesh with constant coloring and saves the video
%   displaced mesh with constant coloring and saves some snapshots
%
% Last modified: February 2, 2018

%Surface refining
pt = [];
ind =1;
for i = 0:k
    pt = [pt; ... 
          (0:k-i)', i*ones(k-i+1,1)];%#ok
    list{i+1} = ind+(0:k-i);%#ok
    ind = ind+k+1-i; % index number for beginning of next row
end

pt = pt/k;
pt = [1-pt(:,1)-pt(:,2), pt];

ell = k+1:-1:1; % ell(i) = length(list{i})

Tri =[];
for i = 1:k
    t = list{i}(1:ell(i)-1); % bottom left vertex
    Tri = [Tri, [t;1+t;ell(i) + t]];%#ok
    t(end) = [];
    Tri = [Tri, [t+1; 1+ell(i)+t; ell(i)+t]];%#ok
end

LT = T.faces(1:3, T.faces(4,:)~=0);
X = T.coordinates(1,:); X=pt*X(LT);
Y = T.coordinates(2,:); Y=pt*Y(LT);
Z = T.coordinates(3,:); Z=pt*Z(LT);

d2 = nchoosek(k+2,2);
TRI = bsxfun(@plus, d2*(0:size(LT,2)-1), Tri(:));
TRI = reshape(TRI, 3, k^2*size(LT,2));

faceList = find(T.faces(4,:)~=0);
DOF=computeBDDOF3D(T,k,faceList);
Nf=size(faceList,2);
% in case of snapshots
iwantsnapshots=false;
if nargin>=5
    snapshots=varargin{2};
    if isnumeric(snapshots)
        uhtot=uhtot(:,snapshots);
    end
    iwantsnapshots=true;
end
P=bernstein2D(pt(:,2),pt(:,3),k);
len = size(uhtot,1)/3;
uh = uhtot(1:len,:); Ntime = size(uh,2); uh=uh(:);
vh = uhtot(len+1:2*len,:); vh=vh(:);
wh = uhtot(2*len+1:end,:); wh=wh(:);
M=ones(d2,1)*(0:Ntime-1); M=len*M(:);
DOFM=bsxfun(@plus,repmat(DOF,[Ntime 1]),M);
uh=uh(DOFM); vh=vh(DOFM); wh=wh(DOFM);
uh=mat2cell(uh,d2*ones(1,Ntime),Nf); uh=cell2mat(uh');
vh=mat2cell(vh,d2*ones(1,Ntime),Nf); vh=cell2mat(vh');
wh=mat2cell(wh,d2*ones(1,Ntime),Nf); wh=cell2mat(wh');
U=P*uh; U=reshape(U,[d2*Nf,Ntime]);
V=P*vh; V=reshape(V,[d2*Nf,Ntime]);
W=P*wh; W=reshape(W,[d2*Nf,Ntime]);

% Taking the values of uh which contributed by the vertices and adding the
% mesh coordinates
Uh=bsxfun(@plus,U,X(:));
Vh=bsxfun(@plus,V,Y(:));
Wh=bsxfun(@plus,W,Z(:));

% Setting a fix coloring
disp('Uniform coloring is set.');
C = ones(Nf,Ntime);
maxc = max(max(C)) + eps; %adding a small number to avoid minc=maxc case
minc = min(min(C));

% Determining the movie options
iwantmovie=false;
rot=0;
angle=[-37.5,30];
facecolor='flat';
axismode=true;
linestyle='-';
framerate=30;
colorscale=parula;
if nargin>=4
    movie=varargin{1};
    iwantmovie=true;
    if ~isempty(movie)
        rot=movie{1};
        if length(movie)>1
            angle=movie{2};
            if strcmp(angle,'default')==1
                angle=[-37.5, 30];
            end
            if length(movie)>2
                facecolor=movie{3};
                if strcmp(facecolor,'default')==1
                    facecolor='flat';
                end
                if ischar(facecolor) && ~isempty(strfind(facecolor,'scale'))%#ok
                    colorscale=strrep(facecolor,'scale','');
                    facecolor='flat';
                end
                if length(movie)>3
                    axismode=movie{4};
                    if length(movie)>4
                        linestyle=movie{5};
                        if length(movie)>5
                            framerate=movie{6};
                        end
                    end
                end
            end
        end
    end
end
if iwantsnapshots
    iwantmovie=false;
end

%creating the movie in the designated folder
if iwantmovie
    close all;
    DateAndTime=datetime('now');
    TimeStamp=datestr(datenum(DateAndTime),'yy-mm-dd-HH-MM');
    if exist([pwd '/movies'],'file') ~= 7
        disp('We can not find the /movies folder in current path.');
        FileStr=strcat(TimeStamp,'_FastMovie_k=',...
            num2str(k),'_M=',num2str(Ntime));
        writerObj = VideoWriter(FileStr,'MPEG-4');
        writerObj.Quality = 100;
        writerObj.FrameRate = framerate;
        open(writerObj);
        disp('Therefore, your videos is being saved in:');
        display(pwd);
    else
        FileStr=strcat('movies/',TimeStamp,'_FastMovie_k=',...
            num2str(k),'_M=',num2str(Ntime));
        writerObj = VideoWriter(FileStr,'MPEG-4');
        writerObj.Quality = 100;
        writerObj.FrameRate = framerate;
        open(writerObj);
        disp('Your videos is being saved in:');
        disp([pwd '/movies']);
    end
end

% Animation part
currentcolor=C(:,1);
plotHandle = @(u,v,w,Col) trisurf(TRI', u, v, w, 0.1+0*u,...
    'LineStyle',linestyle,...
    'FaceLighting','flat',...
    'FaceColor',facecolor,...
    'AmbientStrength',0.3,...
    'DiffuseStrength', 0.8,...
    'SpecularStrength', 0.9,...
    'SpecularExponent', 10);
% Creating a full size figure
figure('units','normalized','outerposition',[0 0 1 1]);
shading interp
colormap(colorscale);
figH=plotHandle(Uh(:,1), Vh(:,1), Wh(:,1), currentcolor);
axis equal;
if ~iwantsnapshots
    title(['time step = ' num2str(1) '/' num2str(Ntime)]);
end
plotmaxx=max(max(Uh)); plotminx=min(min(Uh));
plotmaxy=max(max(Vh)); plotminy=min(min(Vh));
plotmaxz=max(max(Wh)); plotminz=min(min(Wh));
xlim([plotminx,plotmaxx])
ylim([plotminy,plotmaxy])
zlim([plotminz,plotmaxz])
caxis([minc,maxc]);
current_angle = [(1-1)/(Ntime)*rot*360,0] + angle;
%if no rotation is given, angle of view is fixed
%if rotation = 1, then it rotates as the movie goes on
view(current_angle);
if axismode==false
    axis off;
else
    ax = gca;
    ax.XLabel.String = 'x';
    ax.XLabel.FontSize = 12;
    ax.YLabel.String = 'y';
    ax.YLabel.FontSize = 12;
    ax.ZLabel.String = 'z';
    ax.ZLabel.FontSize = 12;
    ax.ZLabel.Rotation = 0;
end
camlight('headlight');

if iwantmovie
    Myvideo.F(1) =  getframe(gcf);
    writeVideo(writerObj,Myvideo.F(1));
elseif iwantsnapshots
    if exist([pwd '/movies'],'file') ~= 7
        figurename=['snapshot_' num2str(1) '.png'];
    else
        figurename=['movies/snapshot_' num2str(1) '.png'];
    end
    saveas(gcf,figurename);
else
    pause(0.1);
end
for j=2:Ntime
    if ~ishghandle(figH)
        break
    end
    shading interp
    colormap(colorscale);
    currentcolor=C(:,j);
    figH=plotHandle(Uh(:,j), Vh(:,j), Wh(:,j), currentcolor);
    axis equal;
    xlim([plotminx,plotmaxx])
    ylim([plotminy,plotmaxy])
    zlim([plotminz,plotmaxz])
    caxis([minc,maxc]);
    if ~iwantsnapshots
        title(['time step = ' num2str(j) '/' num2str(Ntime)]);
    end
    current_angle = [(j-1)/(Ntime)*rot*360,0] + angle;
    %if no rotation is given, angle of view is fixed
    %if rotation = 1, then it rotates as the movie goes on
    view(current_angle);
    if axismode==false
        axis off;
    else
        ax = gca;
        ax.XLabel.String = 'x';
        ax.XLabel.FontSize = 12;
        ax.YLabel.String = 'y';
        ax.YLabel.FontSize = 12;
        ax.ZLabel.String = 'z';
        ax.ZLabel.FontSize = 12;
        ax.ZLabel.Rotation = 0;
    end
    camlight('headlight');
    if iwantmovie
        Myvideo.F(j) =  getframe(gcf);
        writeVideo(writerObj,Myvideo.F(j));
    elseif iwantsnapshots
        if exist([pwd '/movies'],'file') ~= 7
            figurename=['snapshot_' num2str(j) '.png'];
        else
            figurename=['movies/snapshot_' num2str(j) '.png'];
        end
        saveas(gcf,figurename);
    else
        pause(0.1);
    end
end
if iwantmovie
    close(writerObj)
    close
elseif iwantsnapshots
    clc;close all;
    if exist([pwd '/movies'],'file') ~= 7
        disp('Your snapshots are being saved in:');
        disp(pwd);
    else
        disp('Your snapshots are being saved in:');
        disp([pwd '/movies']);
    end
end


end