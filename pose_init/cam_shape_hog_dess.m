function [hogRender] = cam_shape_hog_dess(Shapes, Cameras, Para, verbose)

numShapes = length(Shapes);
numCameras = length(Cameras);
numRenders = numShapes*numCameras;
dimHog = Para.gridHog(1)*Para.gridHog(2)*Para.numOrients*4;
imageDims = Para.gridHog*60;

hogRender = single(zeros(dimHog, numRenders));

for shapeId = 1:numShapes
    if verbose == 1
        fprintf('  Shape_%d\n', shapeId);
        fprintf('    [');
    end
    for cameraId = 1:numCameras
        image = cam_render_shape(Shapes{shapeId}, Cameras{cameraId});
        image = cam_crop_image(image);
        patch = imResample(single(image), imageDims)/255;
        H = hog(patch, 60, Para.numOrients);    
        H = reshape(H, [dimHog,1]);
        hogRender(:, (shapeId-1)*numCameras + cameraId) = H;
        % Compute hog descriptor
        if verbose == 1
           fprintf('.');
           if mod(cameraId, 60) == 0
                fprintf('\n       ');
           end
        end
    end
    if verbose == 1
        fprintf(']\n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that is used to crop images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img_out] = cam_crop_image(img)

[m,n,k] = size(img);
if k > 1
    [Gmag, Gdir] = imgradient(rgb2gray(img), 'roberts');
else
    [Gmag, Gdir] = imgradient(img, 'roberts');
end
Gmag = abs(Gmag);
if m > 1000 || n > 1000
    Gmag(1:4,:) = 0;
    Gmag((m-3):m,:) = 0;
    Gmag(:, 1:4) = 0;
    Gmag(:, (n-3):n) = 0;
end
val = max(max(Gmag));
cols = max(Gmag);
cols = find(cols > val/20);
rows = max(Gmag');
rows = find(rows > val/20);
img = img(rows, cols,:);
img_out = img;

