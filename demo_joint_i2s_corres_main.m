function [Corres_rr, Corres_rs] = demo_joint_i2s_corres_main(...
    Images, ImageCameras, Shapes, Adj_rr, Adj_rs, Para)
% Thie function establishes dense correspondences among a collection of
% real images and natural images
% Input arguments:
%    Images: input images, we assume images are of the same size, and
%            please scale them properly
%    ImageCameras: the associated camera configuration
%    Shapes: input shapes
%    Adj_rr: the sparse matrix that specifies the pairs of images that we
%            want to compute dense correspondences
%    Adj_rs: the sparse matrix that specifies the pairs of image-shapes
%            that we want to compute dense correspondences
%    Para:   please refer to software document on how to change default
%            parameters

% Output arguments:
%    Corres_rr: the dense correspondences between pairs of images
%    Corres_rs: the dense correspondences between images and shapes
%
Clusters = jis_view_clustering(ImageCameras, Para.directDis);
%
for cId = 1:length(Clusters)
    camera = Clusters{cId}.camera;
    imageIds = Clusters{cId}.imageIds;
    %
    numShapes = length(Shapes);
    for id = 1 : numShapes
        [meshPoints, renderImage] = ray_tracing(Shapes{id}, camera);
        images{id} = renderImage;
        scans{id} = meshPoints;
    end
    
end

function [meshPoints, renderImage] = ray_tracing(Shape, Camera)

renderImage = cam_render_shape(Shape, Camera);
[Height, Width, k] = size(renderImage);
cols = kron(1:Width, ones(1, Height));
rows = kron(ones(1, Width), 1:Height);

coordX = ((2*cols-1) - Width)/min(Height, Width);
coordY = (Height - (2*rows-1))/min(Height, Width);
points = axis_x*coordX'*Camera.scale + axis_y*coordY'*Camera.scale +...
    Camera.lookAt*ones(1, numEdgePoints);

meshPoints = unproject(double(Shape.vertexPoss), double(Shape.faceVIds),...
    Camera.origin, points);