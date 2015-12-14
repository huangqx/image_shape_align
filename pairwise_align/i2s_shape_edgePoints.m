function [meshPoints, imagePoints, render_image] = i2s_shape_edgePoints(...
    Shape, Camera, Para)
% Render the shape from a given camera configuration and return indices of
% samples and the rendered image contour
% Input arguments:
%   Shape:  input shape
%   Camera: the camera configuration
% Output argument:
%   meshPoints: coordinates on the original mesh, the last row stores the
%               2D saliency score
%   render_image: rendered color image

% Rectify the camera
axis_z = Camera.origin  - Camera.lookAt;
axis_z = axis_z/norm(axis_z);
Camera.upVec = Camera.upVec - axis_z*(axis_z'*Camera.upVec);
axis_y = Camera.upVec;
axis_x = cross(axis_y, axis_z);

% Rendered 2D image
render_image = i2s_render_shape(Shape, Camera);
edgeMap = edgesDetect(render_image, Para.model_contour);
edgeMap(find(edgeMap < 0.25)) = 0;
edgeMap(find(edgeMap > 0.25)) = 1;
[rows, cols, vals] = find(edgeMap);
numEdgePoints = length(rows);

% Perform ray tracing to find the 3D mesh points
Height = Camera.nHeight_inner;
Width = Camera.nWidth_inner;

coordX = ((2*cols-1) - Width)/min(Height, Width);
coordY = (Height - (2*rows-1))/min(Height, Width);


points = axis_x*coordX'*Camera.scale + axis_y*coordY'*Camera.scale +...
    Camera.lookAt*ones(1, numEdgePoints);

meshPoints = unproject(double(Shape.vertexPoss), double(Shape.faceVIds),...
    Camera.origin, points);

% Perform ray triangle intersection to obtain the 3d locations
ids = find(meshPoints(4,:) < 1e5);
meshPoints = [meshPoints(:, ids); vals(ids)'];
imagePoints = double([coordX(ids)'; coordY(ids)'; vals(ids)']);

