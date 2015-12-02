function [Match] = real_2_render_dense_corres(Image, Shape,...
    Camera)
% This dense correspondences between a real image and rendered images
% Note that in the siggraph paper, this is determined among a collection of
% images and shapes jointly. Please refer to the image_shape_corres package
% for details on a more advanced version
% Input arguments: 
%   Image: The input image object
%   Shape: The shape to be matched
%   Camera: The camera configuration for rendering images

% Render the image 
image = i2s_render_shape(Shape, Camera);
% Compute sift flows


% Compute sift-flow


% Rectify the camera
axis_z = Camera.origin  - Camera.lookAt;
axis_z = axis_z/norm(axis_z);
Camera.upVec = Camera.upVec - axis_z*(axis_z'*Camera.upVec);
axis_y = Camera.upVec;
axis_x = cross(axis_y, axis_z);

Height = Camera.nHeight_inner;
Width = Camera.nWidth_inner;

cols = kron(1:Width, ones(1,Height));
rows = kron(ones(1, Width), 1:Height);

coordX = ((2*cols-1) - Width)/min(Height, Width);
coordY = (Height - (2*rows-1))/min(Height, Width);


points = axis_x*coordX*Camera.scale + axis_y*coordY*Camera.scale +...
    Camera.lookAt*ones(1, Height*Width);

facePartIds = zeros(1, size(Shape.faceVIds, 2));
off = 0;
for i = 1:length(Shape.meshes)
    mesh = Shape.meshes{i};
    nf = size(mesh.faceVIds, 2);
    facePartIds((off+1):(off+nf)) = i;
    off = off + nf;
end

meshPoints = unproject(double(Shape.vertexPoss), double(Shape.faceVIds),...
    Camera.origin, points);
pixelIds = find(meshPoints(4,:) < 10);

pixelFaceIds = zeros(Height, Width);
pixelFaceIds(pixelIds) = facePartIds(meshPoints(1, pixelIds));
