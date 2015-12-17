function [Corres_rs] = demo_joint_i2s_corres_main(...
    Images, ImageCameras, Shapes, Adj_rs, Para)
% Thie function establishes dense correspondences among a collection of
% real images and natural images
% Input arguments:
%    Images: input images, we assume images are of the same size, and
%            please scale them properly
%    ImageCameras: the associated camera configuration
%    Shapes: input shapes
%    Adj_rs: the sparse matrix that specifies the pairs of image-shapes
%            that we want to compute dense correspondences
%    Para:   please refer to software document on how to change default
%            parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output arguments:
%    Corres_rr: the dense correspondences between pairs of images
%    Corres_rs: the dense correspondences between images and shapes
%
Clusters = jis_view_clustering(ImageCameras,...
    Para.viewDirectRadius_image_clustering);
%
fprintf('There are %d clusters.\n', length(Clusters));
%
for cId = 1:length(Clusters)
    fprintf(' processing cluster %d...\n', cId);
    Camera = Clusters{cId}.Camera;
    imageIds = Clusters{cId}.imageIds;
    %
    numShapes = length(Shapes);

    for id = 1 : length(imageIds)
        Clusters{cId}.images{id} = Images{imageIds(id)}.im;
    end
    for id = 1 : numShapes
        [meshPoints, renderImage] = co_ray_tracing(Shapes{id}, Camera);
        Clusters{cId}.images{length(imageIds) + id} = renderImage;
        Clusters{cId}.scans{id} = meshPoints;
    end
    % Compute correspondences among the rendered and the natural images
    [Clusters{cId}.ffds, Clusters{cId}.hogDess] =...
        dense_multi_image_corres(Clusters{cId}.images, Para);
end

numImages = length(Images);
imageClusIds = zeros(2, numImages);
for cId = 1 : length(Clusters)
    imageClusIds(1, Clusters{cId}.imageIds) = cId;
    imageClusIds(2, Clusters{cId}.imageIds) =...
        1:length(Clusters{cId}.imageIds);
end

[imIds, shapeIds, vals] = find(Adj_rs);

Corres_rs = cell(1, length(imIds));
for id = 1:length(imIds)
    Corres_rs{id}.imageId = imIds(id);
    Corres_rs{id}.shapeId = shapeIds(id);
    cId = imageClusIds(1, imIds(id));
    im_cId = imageClusIds(2, imIds(id));
    Corres_rs{id}.corres = establish_2d_3d_corres(...
        Clusters{cId}.ffds{im_cId},...
        Clusters{cId}.ffds{shapeIds(id)+length(Clusters{cId}.imageIds)},...
        Clusters{cId}.scans{shapeIds(id)});
end


