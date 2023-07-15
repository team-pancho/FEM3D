function meshPlot(T, f, varargin )
%meshPlot(T, f )
%meshPlot(T, f, 1)   
% Inputs:
%       T : an enhanced tetrahedrization 
%       f : either a function of three variables or a row vector which 
%           is used as the colormap of the plot
%       varargin : type 1 to plot the normals
% 
% Last Modified: November 4, 2016

if size(f,2) > 1
    tetramesh(T.elements',T.coordinates',f);
else
    bary = (T.coordinates(:, T.elements(1,:))...
          + T.coordinates(:, T.elements(2,:))...
          + T.coordinates(:, T.elements(3,:))...
          + T.coordinates(:, T.elements(4,:)))/4; 
    tetramesh(T.elements',T.coordinates',f(bary(1,:),bary(2,:),bary(3,:)));
end

if nargin == 3 && varargin{1} == 1
    baryFace = (T.coordinates(:,T.faces(1,:))...
              + T.coordinates(:,T.faces(2,:))...
              + T.coordinates(:,T.faces(3,:)))/3;
    hold on
    quiver3(baryFace(1,:),baryFace(2,:),baryFace(3,:),...
            T.normals(1,:),T.normals(2,:),T.normals(3,:));
    hold off
end
    

end

