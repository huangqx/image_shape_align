function [renderImage, pixelPartIds, meshPoints] = mm_render_shape(Shape, Camera)
% Assuming that the Shape consists of a set of parts, this function render
% the shape from the given camera configuration 
% Render the shape using the camera configuration
renderImage = cam_render_shape(Shape, Camera);
% Compute sift flows

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
pixelIds = find(meshPoints(4,:) < 1e5);

pixelPartIds = zeros(Height, Width);
pixelPartIds(pixelIds) = facePartIds(meshPoints(1, pixelIds));
