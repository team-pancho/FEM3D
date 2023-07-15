function plotElasticity3D(uhtot,T,k,varargin)
% plotElasticity3D(uhtot,T,k)
% plotElasticity3D(uhtot,T,k,sigmah)
% plotElasticity3D(uhtot,T,k,sigmah,{})
%     OR     plotElasticity3D(uhtot,T,k,sigmah,movie)
% plotElasticity3D(uhtot,T,k,sigmah,movie,snapshots)
%
% Input:
%   uhtot    :      3*dimVh x Ntimes complex vector (containing uh,vh,wh)
%   T        :      Data structure. Enhanced FEM triangulation
%   k        :      Scalar. Polynomial degree
%   sigmah   :      Nelts x Ntimes numerical norm of stress at barycenter
%                   of elements. If stress is unknown, and uniform coloring
%                   is desired pass in a scalar i.e. 1.
%   movie    :      Cell element containing the following
%                   (for default you can pass in {})
%   movie{1} :      Rotation. This can be 1 or 0, enables or disables
%                     camera rotation.
%                     Default: no rotation
%   movie{2} :      Camera angle. 1x2 vector corresponding to
%                     [azimuth elevation] in radians.
%                     Input option: 'default'
%                     Default: [-37.5,30]
%   movie{3} :      Color for the faces in MATLAB RGB vector,
%                     Input options: 'default' OR 'flat',
%                      OR 'grayscale', 'hotscale, ''winterscale',...
%                      OR 'green','red',...
%                      OR [1,0.3,0.8]
%                     'flat': Uniform colors taken from stress data
%                     'grayscale': Computes the grayscale of the
%                                  given coloring in the stress data
%                     'green': Uniform red color ignoring stress
%                     Default: 'flat'
%   movie{4} :      Axis on or off. This can be 1 or 0. When
%                     it is on it shows the labels too.
%                     Default: on
%   movie{5} :      Edge line style.
%                     Input options: '-'(solid line), 'none', '--'(dashed),
%                     ':' (dotted)
%                     Default: '-'
%   movie{6} :      Number of frames per second. A positive integer.
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
%
% Action: plot of ...
%   displaced mesh with constant coloring
%   displaced mesh colored by norm of stress
%   displaced mesh colored by norm of stress and saves the video
%   displaced mesh colored by norm of stress and saves the snapshots
%
% Last modified: February 1, 2018

% Finding vertices from DOF
Ne = size(T.elements,2);
dofuh = DOF3D(T,k);
vdof = unique(dofuh(1:4,:));
% in case of snapshots
iwantsnapshots=false;
if nargin>=6
    snapshots=varargin{3};
    if isnumeric(snapshots)
        uhtot=uhtot(:,snapshots);
    else
        snapshots=1:size(uhtot,2);
    end
    iwantsnapshots=true;
end
% Decomposing uhtot
len = size(uhtot,1)/3;
uh = uhtot(1:len,:);
vh = uhtot(len+1:2*len,:);
wh = uhtot(2*len+1:end,:);
Ntime = size(uh,2);
% Taking the values of uh which contributed by the vertices and adding the
% mesh coordinates
Uh=bsxfun(@plus,uh(vdof,:),T.coordinates(1,:)');
Vh=bsxfun(@plus,vh(vdof,:),T.coordinates(2,:)');
Wh=bsxfun(@plus,wh(vdof,:),T.coordinates(3,:)');

% Finding face numbering
faceList=[T.dirichlet T.neumann];
Nf=size(faceList,2);
[~,~,ib]=intersect(faceList.',T.faces(1:3,:).','rows','stable');
[~,~,ibb]=intersect(ib,T.facebyelt(:),'stable');
colno=(ceil(ibb/4));
% Determining the coloring via using stress
if nargin>=4
    C = varargin{1};
    if size(C,1)==Ne
        C=C(colno,:);
        if nargin>=6
            C=C(:,snapshots);
        end
        disp('Face coloring is derived from given col set.');
    end
    if size(C,2)~=Ntime
        disp('Uniform coloring is set.');
        C = ones(Nf,Ntime);
        if nargin>6
            C=C(:,snapshots);
        end
    end
else
    disp('Uniform coloring is set.');
    C = ones(Nf,Ntime);
end
maxc = max(max(C)) + eps; %adding a small number to avoid minc=maxc case
minc = min(min(C));

% Determining the movie options
iwantmovie=false;
rot=0;
angle=[-37.5,30];
facecol='flat';
axismode=true;
linestyle='-';
framerate=30;
colorscale=parula;
if nargin>=5
    movie=varargin{2};
    iwantmovie=true;
    if ~isempty(movie)
        rot=movie{1};
        if length(movie)>1
            angle=movie{2};
            if strcmp(angle,'default')==1
                angle=[-37.5, 30];
            end
            if length(movie)>2
                facecol=movie{3};
                if strcmp(facecol,'default')==1
                    facecol='flat';
                end
                if ischar(facecol) && ~isempty(strfind(facecol,'scale'))%#ok
                    colorscale=strrep(facecol,'scale','');
                    facecol='flat';
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
plotHandle = @(u,v,w,Col) trisurf(faceList', u, v, w, Col,...
    'LineStyle',linestyle,...
    'FaceLighting','flat',...
    'FaceColor',facecol,...
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
        return
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