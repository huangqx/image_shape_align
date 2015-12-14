function [Shape_opt, Camera_opt] = image2shape_part_based_align(...
    Image,... % The input image object
    Camera_init,... % The initial camera configuration
    Shape_init,... % The initial assembled shape
    partCorres,... % The part correspondences (currently un-used)
    Para)
% This function aligns the input shape with the image object
% Input arguments:
%   Image: The input image object
%   Camera_init: The initial camera configuration
%   Shape_init: The initial shape assembed from a set of parts
%   partCorres: The correspondences between image pixels and part mesh
%               points
%   Para: Please refer to the software document
% Output arguments:
%   Shape_opt: The optimized shape
%   Camera_opt: The associated optimized camera configuration

% Perform part-based fine-tunning
% Compute image edges                % set to true to enable nms
edgeMap = edgesDetect(Image.im, Para.model_contour);

% Crop the image based on bounding boxes
edgeMap(:, 1:floor(Image.bbox(1))) = 0;
edgeMap(:, floor(Image.bbox(3)):size(Image.im,2)) = 0;
edgeMap(1:floor(Image.bbox(2)),:) = 0;
edgeMap(floor(Image.bbox(4)):size(Image.im,1),:) = 0;

% Crop non-salient pixels
if 1
    edgeMap(find(edgeMap < 0.25)) = 0;
    edgeMap(find(edgeMap > 0.25)) = 1;
end

% Extract boundary pixels
[rows, cols, vals] = find(edgeMap);
numPixels = length(rows);
targetPC_2D = zeros(3, numPixels);
[Height, Width] = size(edgeMap);

% Image contours to be aligned
targetPC_2D(1,:) = ((2*cols-1) - Width)/min(Height, Width);
targetPC_2D(2,:) = (Height - (2*rows-1))/min(Height, Width);
targetPC_2D(3,:) = vals;


for i = 1:length(partCorres)
    partCorres{i}.points2d(1,:) = ((2*partCorres{i}.pixels(2,:)-1) - Width)/min(Height, Width);
    partCorres{i}.points2d(2,:) = (Height - (2*partCorres{i}.pixels(1,:)-1))/min(Height, Width);
    %
    mesh = Shape_init.meshes{i};
    mesh.vertexPoss = Shape_init.vertexPoss(:, mesh.vertexIds);
    %
    partCorres{i}.points3d = extact_3d_points(mesh,...
        partCorres{i}.meshPoints);
    partCorres{i}.weights = ones(1, size(partCorres{i}.points3d, 2));
end


Shape_opt = Shape_init;
Camera_opt = Camera_init;

partTrans = alloc_part_trans(Shape_init);

% Perform ICP style optimization
for iter = 1:Para.numIterations
    fprintf(' iteration_%d\n', iter);
    % Compute feature points of a rendered image and locate them in 3D
    fprintf('   building correspondences\n');
    [meshPoints, sourcePC_2D, render_image] = i2s_shape_edgePoints(...
        Shape_opt,...
        Camera_opt,...
        Para);
    
    if 0 % Used for debugging
        hFig = figure(1);
        imshow((im2double(Image.im) + im2double(render_image))/2);
        close(hFig);
    end
    
    
    
    
    % Given the mesh coordinates, find the 3D locations
    if iter > 2
        Corres = i2s_nn_query_pruned(...
        targetPC_2D,...
        sourcePC_2D,...
        Para.minSigmaCorres,...
        Para.lambda_saliency);
    
        partCorres = extract_part_corres(Shape_init,...
            targetPC_2D,...
            meshPoints,...
            Corres);
    end
    
    fprintf('   part-based alignment...\n');
    % Current part-based FFD is deactivated
    [Camera_opt, partTrans] = mm_part_affine_align(...
        Shape_init,...
        partTrans,...
        partCorres,...
        Camera_opt,...
        Para);
    
    % Transform the shape
    Shape_opt = transform_shape(Shape_init, partTrans);
end


if 1 % Used for debugging
    hFig = figure(2);
    imshow((im2double(Image.im) + 2*im2double(render_image))/3);
end

function [Shape_opt] = transform_shape(Shape, partTrans)
%
Shape_opt = Shape;
numSyms = length(Shape.nonsym_ids);
numPairs = size(Shape.sym_pairs, 2);
for i = 1:numSyms
    partId = Shape.nonsym_ids(i);
    theta = partTrans.theta_sym(i);
    R = [1, 0, 0;
        0, cos(theta), -sin(theta);
        0, sin(theta), cos(theta)];
    R = eye(3);
    t = [0, partTrans.Trans_yz(1, i), partTrans.Trans_yz(2,i)]';
    vertexPoss = Shape.vertexPoss(:, Shape.meshes{partId}.vertexIds);
    Shape_opt.vertexPoss(:, Shape_opt.meshes{partId}.vertexIds) =...
        R*vertexPoss + t*ones(1, size(vertexPoss,2));
end
for i = 1:numPairs
    part1Id = Shape.sym_pairs(1, i);
    R = partTrans.Rot_pairs(:,:,i);
    R = eye(3);
    t = partTrans.Trans_pairs(:,i);
    vertexPoss = Shape.vertexPoss(:,...
        Shape.meshes{part1Id}.vertexIds);
    Shape_opt.vertexPoss(:, Shape_opt.meshes{part1Id}.vertexIds) =...
        R*vertexPoss + t*ones(1, size(vertexPoss,2));
 
    part2Id = Shape.sym_pairs(2, i);
    vertexPoss = Shape.vertexPoss(:,...
        Shape.meshes{part2Id}.vertexIds);
        
    t(1) = -t(1);
    R = diag([-1,1,1])*R*diag([-1,1,1]);
    Shape_opt.vertexPoss(:, Shape_opt.meshes{part2Id}.vertexIds) =...
        R*vertexPoss + t*ones(1, size(vertexPoss,2));
end


function [Var] = alloc_part_trans(Shape)
%
numPairs = size(Shape.sym_pairs, 2);
numSyms = length(Shape.nonsym_ids);
Var.Rot_pairs = zeros(3, 3, numPairs);
Var.Trans_pairs = zeros(3, numPairs);
for i = 1:numPairs
    Var.Rot_pairs(:,:,i) = eye(3);
    Var.Trans_pairs(:,i) = zeros(3, 1);
end
Var.theta_sym = zeros(1, numSyms);
Var.Trans_yz = zeros(2, numSyms);


function [partCorres] = extract_part_corres(Shape_ori,...
    targetPC_2D,...
    meshpoints,...
    Corres)
% extract part-level correspondences
nf = size(Shape_ori.faceVIds, 2);
facePartIds = zeros(1, nf);
nf = 0;
for partId = 1:length(Shape_ori.meshes)
    mesh = Shape_ori.meshes{partId};
    nf_mesh = size(mesh.faceVIds, 2);
    facePartIds((nf+1):(nf+nf_mesh)) = partId;
    nf = nf + nf_mesh;
end

points2d = targetPC_2D(:, Corres(2,:));
tp = extact_3d_points(Shape_ori, meshpoints);
points3d = tp(:, Corres(1,:));
meshpoints = meshpoints(:, Corres(1,:));
weights = Corres(3, :);

facePartIds = facePartIds(meshpoints(1,:));

numParts = length(Shape_ori.meshes);
partCorres = cell(1, numParts);
for partId = 1 : numParts
    ids = find(facePartIds == partId);
    partCorres{partId}.points2d = points2d(:, ids);
    partCorres{partId}.points3d = points3d(:, ids);
    partCorres{partId}.weights = weights(ids);
end

function [points3d] = extact_3d_points(Shape, meshpoints)


faceVIds = Shape.faceVIds(:, meshpoints(1,:));
points_1 = Shape.vertexPoss(:, faceVIds(1,:));
points_2 = Shape.vertexPoss(:, faceVIds(2,:));
points_3 = Shape.vertexPoss(:, faceVIds(3,:));

paras1 = 1 - meshpoints(2,:) - meshpoints(3,:);
paras2 = meshpoints(2,:);
paras3 = meshpoints(3,:);
points3d = points_1.*(ones(3,1)*paras1)+...
    points_2.*(ones(3,1)*paras2)+...
    points_3.*(ones(3,1)*paras3);
