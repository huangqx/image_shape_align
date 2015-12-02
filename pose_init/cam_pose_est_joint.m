function [optimizedPoses] = cam_pose_est_joint(...
    inputImages,...
    Shapes,...
    cameraSamples,...
    hogRender,...
    Para,...
    verbose)
% Estimate the camera pose of a collection of images given a shape
% collection
% Input argument:
%       inputImages: the input image objects
%       Shapes: a cell array that stores all the shape objects
%       cameraSamples: the discrete set of camera samples, the camera
%               pose of each image will be selected from this array
%       Para: stores all the parameters used in this step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters used in computing the HOG descriptor
% Para.gridHog = [6, 6]
% Para.numOrients = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numImages = length(inputImages);
numShapes = length(Shapes);
numCameras = length(cameraSamples);
numRenders = numShapes*numCameras;
if numRenders ~= size(hogRender, 2)
    fprintf('The dimension of the rendered images has to agree with\n');
    fprintf(' the dimension of the descriptor.\n');
end
dimHog = Para.gridHog(1)*Para.gridHog(2)*Para.numOrients*4;
imageDims = Para.gridHog*60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute hog descriptors for the input images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose == 1
    fprintf('Compute desctiptors for input images...\n');
end

hogInput = single(zeros(dimHog, numImages));
for imId = 1 : numImages
    bbox = floor(inputImages{imId}.bbox);
    image = inputImages{imId}.im(bbox(2):bbox(4),...
        bbox(1):bbox(3),:);
    patch = imResample(single(image), imageDims)/255;
    H = hog(patch, 60, Para.numOrients);    
    H = reshape(H, [dimHog,1]);
    hogInput(:, imId) = H;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the unary term for the MRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unaryTerm = zeros(numCameras, numImages);

for imId = 1 : numImages
    hog_im = hogInput(:, imId);
    dif = hog_im*single(ones(1, numRenders)) - hogRender;
    dif = reshape(sqrt(sum(dif.*dif)), [numCameras, numShpaes]);
    dif = sort(dif')';
    scores = dif(:, 1:Para.knn_unary);
    unaryTerm(:, imId) = mean(scores')';
end

fprintf('Compute image-graph\n');
W_cam = cam_pose_similarity_matrix(cameraSamples);
simGraph = cam_image_graph(hogInput, knn);
% Joint inference of the camera poses
% 
sol = cam_joint_inference(unaryTerm,...
    simGraph,...
    W_cam,...
    Para.lambda_graph);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Perform pose estimation\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for imId = 1 : numImages
    cameraId = sol(imId);
    optimizedPoses{imId} = cameraSamples{cameraId};
    hog_im = hogInput(:, imId);
    dif = hog_im*single(ones(1, numRenders)) - hogRender;
    dif = reshape(sqrt(sum(dif.*dif)), [numCameras, numShpaes]);
    dif = dif(cameraId, :);
    [scores, shapeIds] = sort(dif);
    optimizedPoses{imId}.closestShapeIds =...
        shapeIds(1:Para.knn_unary);
end

function [W_cam] = cam_pose_similarity_matrix(cameraSamples)
%
numCameras = length(cameraSamples);
poss = zeros(3, numCameras);
for id = 1:numCameras
    poss(:, id) = Cameras{id}.origin;
end

nnDis = zeros(numCameras, 1);
for id = 1:numCameras
    dif = poss(:, id)*ones(1, numCameras) - poss;
    dif = sqrt(sum(dif.*dif));
    [s,id] = sort(dif);
    nnDis(id) = s(2);
end
sigma = mean(nnDis);

rowIds = kron(1:numCameras, ones(1,numCameras));
colIds = kron(ones(1,numCameras), 1:numCameras);
dif = poss(:, rowIds) - poss(:, colIds);
dif = sqrt(sum(dif.*dif));
dif = exp(-dif.*dif/2/sigma/sigma);
W_cam = reshape(dif, [numCameras, numCameras]);
%

function [sol] = cam_joint_inference(unaryTerm,...
    imageGraph,...
    W_cam,...
    lambda_graph)
% write the interface to opengm2
%
[numStates, numImages] = size(unaryTerm);
offs = 0:numStates:(numStates*numImages);
stateDims = ones(1, numImages)*numStates;
%
unaryData = reshape(unaryTerm, [1, numStates*numImages]);
%
[sVIds, tVIds, edgeWeights] = find(imageGraph);
numEdges = length(sVIds);
%
dim2 = numStates*numStates;
w_cam = reshape(W_cam, [1, dim2]);
binaryData = kron(ones(1, numEdges), w_cam);
binaryData = lambda_graph*binaryData.*kron(ones(1,dim2), edgeWeights');

paras = [400, 1e-9, 0.25];
sol = trws(stateDims, unaryData, sVIds'-1, tVIds'-1, -binaryData, paras);


function [G] = cam_image_graph(hogInput, knn)

numImages = size(hogInput, 2);
rowsG = (1:numImages)'*ones(1, knn);
colsG = zeros(numImages, knn);
valsG = zeros(numImages, knn);

G = sparse(numImages, numImages);
for sId = 1:numImages
    dif = hogInput(:, sId)*single(ones(1, numImages)) - hogInput;
    dif = sqrt(sum(dif.*dif));
    [s, ids] = sort(dif);
    colsG(sId,:) = s(1:knn);
    valsG(sId,:) = ids(1:knn);
end

rowsG = reshape(rowsG, [1, knn*numImages]);
colsG = reshape(colsG, [1, knn*numImages]);
valsG = reshape(valsG, [1, knn*numImages]);

sigma = median(valsG);
G = sparse(rowsG, colsG, exp(-(valsG.*valsG)/2/sigma/sigma),...
    numImages, numImages);
G = max(G, G');