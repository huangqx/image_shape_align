function [knnIds] = cam_knn_nearest_query(Image, Camera, Shapes, knn)
% This function finds the nearest neighbors of an image object among a
% collection of shapes stored in a shape folder. The shapes are assumed to
% be consistently oriented and scaled.
% Input arguments:
%   Image: the input image object
%   Camera: the associated camera configuration
%   shape_folder: the shape folder that stores all the input shapes, The
%                 shapes are assumed to be consistently oriented, and the
%                 length of the diagonal of each bounding box is 1 
%                 (the default setting in pose estimation) 
%   knn: the number of nearest-neighbors for a given image object
% Output argument:
%   knnIds: indices of the retrieved nearest neighbors

patch = imResample(single(Image.im), [360, 360])/255;
H = hog(patch, 60, 16);
query_hog = reshape(H, [2304, 1]);

numShapes = length(Shapes);

shapes_hog = single(zeros(2304, numShapes));


for shapeId = 1 : numShapes
    shapes_hog(:, shapeId) = hos_process(Shapes{shapeId}, Camera);
    if mod(shapeId, 100) == 0
        fprintf('%d\n', shapeId);
    end
end


d = query_hog*ones(1, numShapes) - shapes_hog;
[s, ids] = sort(sum(d.*d));
knnIds = ids(1:knn);


function [H] = hos_process(Shape, Camera)

Shape.has_material = 1;
image = cam_render_shape(Shape, Camera);
patch = imResample(single(image), [360, 360])/255;
H = hog(patch, 60, 16);
H = reshape(H, [2304, 1]);