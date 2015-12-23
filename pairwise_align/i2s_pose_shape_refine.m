function [Camera_opt, Shape_opt, FFD_opt] = i2s_pose_shape_refine(...
    Image, Shape_in, Camera_in, Para)
% Image: the pascal3d image to be aligned
% Shape: the template surface

% Initialize the normalize the cameras
Camera_opt = Camera_in;
Camera_opt = i2s_set_camera(Image.im, Camera_opt,...
    Para.rHeights, Para.rWidths);
Camera_opt = i2s_normalize_camera(Image.bbox,...
    Shape_in, Camera_opt);

% Compute image edges                % set to true to enable nms
edgeMap = edgesDetect(Image.im, Para.model_contour);

% Crop the image based on bounding boxes
edgeMap(:, 1:floor(Image.bbox(1))) = 0;
edgeMap(:, floor(Image.bbox(3)):size(Image.im,2)) = 0;
edgeMap(1:floor(Image.bbox(2)),:) = 0;
edgeMap(floor(Image.bbox(4)):size(Image.im,1),:) = 0;

% Crop non-salient pixels
edgeMap(find(edgeMap < 0.15)) = 0;
[rows, cols, vals] = find(edgeMap);

numPixels = length(rows);
targetPC_2D = zeros(3, numPixels);
[Height, Width] = size(edgeMap);

% Image contours to be aligned
targetPC_2D(1,:) = ((2*cols-1) - Width)/min(Height, Width);
targetPC_2D(2,:) = (Height - (2*rows-1))/min(Height, Width);
targetPC_2D(3,:) = vals;

Shape_opt = Shape_in;

% Precomputed ffd coefficients for shape vertices and shape samples
Shape_in.FFD = i2s_ffd_init_sym(Shape_in, Para.ffd_Res);
% ffd_coeff_edgeSam = i2s_ffd_basis_coeff(Shape_in.FFD,...
%     double(Shape_in.Bound.edgeSamples));
% ffd_coeff_vertexSam = i2s_ffd_basis_coeff(Shape_in.FFD,...
%     double(Shape_in.Bound.vertexSamples));
ffd_coeff_vertexPoss = i2s_ffd_basis_coeff(Shape_in.FFD,...
    double(Shape_in.vertexPoss));

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
    
    
    Corres = i2s_nn_query_pruned(...
        targetPC_2D,...
        sourcePC_2D,...
        Para.minSigmaCorres,...
        Para.lambda_saliency);
    
    
    % Given the mesh coordinates, find the 3D locations
    sourcePC_3D_ori = extact_3d_points(Shape_in, meshPoints);
    
    if iter <= Para.numIterations/2
        fprintf('   rigid alignment...\n');
        Camera_opt = i2s_camera_opt(...
            targetPC_2D(:, Corres(2,:)),...
            sourcePC_3D_ori(:, Corres(1,:)),...
            Corres(3,:),...
            Camera_opt,...
            Para);
    else
        fprintf('   non-rigid alignment...\n');
        [Camera_opt, Shape_in.FFD.ctrlPos_cur] =...
            i2s_camera_shape_opt(...
            targetPC_2D(:, Corres(2,:)),...
            sourcePC_3D_ori(:, Corres(1,:)),...
            Corres(3,:),...
            Shape_in.FFD,...
            Camera_opt,...
            Para);
         
        Shape_opt.vertexPoss = single(Shape_in.FFD.ctrlPos_cur*...
            ffd_coeff_vertexPoss');
    end
end

FFD_opt = Shape_in.FFD;

if 1 % Used for debugging
    hFig = figure(2);
    imshow((im2double(Image.im) + 2*im2double(render_image))/3);
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
