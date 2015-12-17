function [meshPoints, renderImage] = co_ray_tracing(Shape, Camera)
%
axis_z = Camera.origin - Camera.lookAt;
viewDis = norm(axis_z);
axis_z = axis_z/viewDis;
axis_y = Camera.upVec;
axis_x = cross(axis_y, axis_z);
%
renderImage = cam_render_shape(Shape, Camera);
[Height, Width, k] = size(renderImage);
cols = kron(1:Width, ones(1, Height));
rows = kron(ones(1, Width), 1:Height);
numPoints = Height*Width;

coordX = ((2*cols-1) - Width)/min(Height, Width);
coordY = (Height - (2*rows-1))/min(Height, Width);
points = axis_x*coordX*Camera.scale + axis_y*coordY*Camera.scale +...
    Camera.lookAt*ones(1, numPoints);

meshPoints = unproject(double(Shape.vertexPoss), double(Shape.faceVIds),...
    Camera.origin, points);
%
mask = rgb2gray(renderImage) < 240;
[rows, cols, vals] = find(mask);
row_min = min(rows);
row_max = max(rows);
col_min = min(cols);
col_max = max(cols);
IDX = reshape(1:(Height*Width), [Height, Width]);
renderImage = renderImage(row_min:row_max, col_min:col_max,:);
IDX = IDX(row_min:row_max, col_min:col_max);
IDX = reshape(IDX, [1, size(IDX,1)*size(IDX,2)]);
meshPoints = meshPoints(:, IDX);