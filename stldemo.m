%% 3D Model Demo
% This is short demo that loads and renders a 3D model of a human femur. It
% showcases some of MATLAB's advanced graphics features, including lighting and
% specular reflectance.

% Copyright 2011 The MathWorks, Inc.


%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
clear all
clc
%fv = stlread('femur.stl');
fv = stlread('tube2_n.stl');
%S=skeleton(fv.vertices);

%% extract centerline z=26'
for i=0:50
    if(i<20)
        x(i+1)=13;
        z(i+1)=13;
        y(i+1)=i;
    else 
        z(i+1)=13;
        y(i+1)=i;
        x(i+1)=-sqrt(30*30-(y(i+1)-20)^2)+43;
    end
end



%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.
%facecolor original [0.8 0.8 1.0]
patch(fv,'FaceColor',       [0 0.9 0.2], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15,           ...
        'FaceAlpha',.2,'EdgeAlpha',.3 );
hold on
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
%camlight('right');
%material('dull');
material('shiny');

plot3(x,y,z,'r','LineWidth',3)
% Fix the axes scaling, and set a nice view angle
axis('image');
xlabel('X(mm)','FontSize',12)
ylabel('Y(mm)','FontSize',12)
zlabel('Z(mm)','FontSize',12)
view([-135 35]);