function [Shapes_knn] = cam_nearest_query(Image, Camera, shape_folder, knn)
% This function finds the nearest neighbors of an image object among a
% collection of shapes stored in a shape folder. The shapes are assumed to
% be consistently oriented and scaled.

patch = imResample(single(Image.im), [360, 360])/255;
H = hog(patch, 60, 16);
query_hog = reshape(H, [2304, 1]);

temp = dir(shape_folder);
numShapes = length(temp)-2;

shapes_hog = single(zeros(2304, numShapes));


for shapeId = 1 : numShapes
    shapes_hog(:, shapeId) = hos_process(shape_folder, temp, Camera, shapeId);
    fprintf('%d\n', shapeId);
end


d = query_hog*ones(1, numShapes) - shapes_hog;
[s, ids] = sort(sum(d.*d));
ids = ids(1:knn);

figure(1);
for i = 1:knn
    shapeId = ids(i);
    load([shape_folder, temp(shapeId+2).name]);
    Shape.has_material = 1;
    Shapes_knn{i} = Shape;
    images{i} = cam_render_shape(Shape, Camera);
    subplot(6,8, i);
    imshow(images{i});
end

function [H] = hos_process(shape_folder, temp, Camera, shapeId)

load([shape_folder, temp(shapeId+2).name]);
Shape.has_material = 1;
image = cam_render_shape(Shape, Camera);
patch = imResample(single(image), [360, 360])/255;
H = hog(patch, 60, 16);
H = reshape(H, [2304, 1]);