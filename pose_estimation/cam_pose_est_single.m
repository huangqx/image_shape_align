function [Camera, im] = cam_pose_est_single(Image, Shapes, Cameras, hogRender,...
    Para, verbose)
% Estimate the camera pose of the input images given the shapes
% Image: the input image object
% Shapes: a cell array that stores all the shape objects
% Para: stores all the parameters used in this step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters associated with the camera sampling
% Parameters used in computing the HOG descriptor
% Para.gridHog = [6, 6]
% Para.numOrients = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sample candidate camera configurations
if verbose == 1
    fprintf('Sampling cameras...\n');
end

if verbose == 1
    fprintf('Render shapes and compute desctiptors...\n');
end

numShapes = length(Shapes);
numCameras = length(Cameras);
numRenders = numShapes*numCameras;
dimHog = Para.gridHog(1)*Para.gridHog(2)*Para.numOrients*4;
imageDims = Para.gridHog*60;

%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute hog descriptors for the input images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose == 1
    fprintf('Compute desctiptors for input images...\n');
end


bbox = floor(Image.bbox);
image = Image.im(bbox(2):bbox(4), bbox(1):bbox(3),:);
patch = imResample(single(image), imageDims)/255;
H = hog(patch, 60, Para.numOrients);    
H = reshape(H, [dimHog,1]);
hogInput = H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the unary term for the MRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To accerlate the process, we do dimension over the
% descriptor space
%projMat = pe_hog_dim_reduct_matrix(hogInput, hogRender);
%hogInput = projMat*hogInput;
%hogRender = projMat*hogRender;

minScore = 1e10;
bestCameraId = 0;
for cameraId = 1:numCameras
    hogThisView = hogRender(:, cameraId:numCameras:numRenders);
    dif = hogInput*single(ones(1, numShapes)) - hogThisView;
    dif = sqrt(sum(dif.*dif));
    [dis, ids] = sort(dif);
    score = mean(dis(1:Para.knn_unary));
    if score < minScore
        minScore = score;
        bestCameraId = cameraId;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick the closest shapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute image-graph\n');
fprintf('Perform pose estimation\n');
figure(1);
hogThisView = hogRender(:, bestCameraId:numCameras:numRenders);
dif = hogInput*single(ones(1, numShapes)) - hogThisView;
Camera = Cameras{bestCameraId};
[s, closetShapeId] = min(sum(dif.*dif));
Camera.closetShapeId = closetShapeId;

% Set the camera
Camera.nHeight_inner = size(Image.im, 1);
Camera.nWidth_inner = size(Image.im, 2);
tp = find(Para.rHeights == Camera.nHeight_inner);
Camera.nHeight = tp(1);
tp = find(Para.rWidths == Camera.nWidth_inner);
Camera.nWidth = tp(1);


im11 = cam_crop_image(Image.im);
im2 = cam_render_shape(Shapes{closetShapeId}, Camera);
im21 = cam_crop_image(im2);
[m1, n1, k] = size(im11);
[m2, n2, k] = size(im21);
s1 = norm([m1, n1]);
s2 = norm([m2, n2]);
Camera.scale = Camera.scale*(s2/s1);
if verbose == 1
    im2 = cam_render_shape(Shapes{closetShapeId}, Camera);
    [im] = stitch_image(Image.im, im2);
    figure(2);
    imshow(im);
end


function [im] = stitch_image(im1, im2)
[m1,n1,k] = size(im1);
[m2,n2,k] = size(im2);
s1 = 500/m1;
s2 = 500/m2;
im1 = imresize(im1, [floor(m1*s1+0.1), floor(n1*s1+0.1)]);
im2 = imresize(im2, [floor(m2*s2+0.1), floor(n2*s2+0.1)]);
n = size(im1, 2) + size(im2, 2);
im = 255 - uint8(zeros(500, n, 3));
im(:,1:size(im1,2),:) = im1;
im(:,(size(im1,2)+1):n,:) = im2;



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

