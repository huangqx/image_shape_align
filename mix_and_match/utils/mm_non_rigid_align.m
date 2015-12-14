function [Shape_opt, Camera_opt] = mm_non_rigid_align(Image, Shape_init,...
    Camera_init, partCorres, Para)
% Align the input shape to the input image object with given
% correspondences. The deforming to the input shape is a free-from
% deformation
% Input arguments:
%      Image: The input image
%      Shape_init: the initial assemble shape
%      Camera_init: the initial camera configuration
%      partCorres: part-wise dense correspondences
% Output arguments:
%      Shape_opt: optimized shape
%      Camera_opt: optimized camera configuration

% Convert partCores into dense point-wise 2D-3D correspondences
numParts = length(partCorres);
numCorres = 0;
for partId = 1:numParts
    pc = partCorres{partId};
    numCorres = numCorres + size(pc.meshPoints, 2);
end
points2d = zeros(2, numCorres);
points3d = zeros(3, numCorres);

numCorres = 0;
for partId = 1:numParts
    pc = partCorres{partId};
    mesh = Shape_init.meshes{partId};
    vertexPoss = Shape_init.vertexPoss(:, mesh.vertexIds);
    p1 = vertexPoss(:, mesh.faceVIds(1, pc.meshPoints(1,:)));
    p2 = vertexPoss(:, mesh.faceVIds(2, pc.meshPoints(1,:)));
    p3 = vertexPoss(:, mesh.faceVIds(3, pc.meshPoints(1,:)));
    
    t1 = ones(3,1)*(1 - pc.meshPoints(2,:) - pc.meshPoints(3,:));
    t2 = ones(3,1)*pc.meshPoints(2, :);
    t3 = ones(3,1)*pc.meshPoints(3, :);
    
    nf = size(pc.meshPoints, 2);
    points2d(:, (numCorres+1):(numCorres+nf)) = pc.pixels;
    points3d(:, (numCorres+1):(numCorres+nf)) = p1.*t1 + p2.*t2 + p3.*t3;
    numCorres = numCorres + nf;
end


corresWeights = ones(1, numCorres);

[Height, Width, k] = size(Image.im);


% Image contours to be aligned
targetPC_2D = zeros(3, numCorres);
targetPC_2D(1,:) = ((2*points2d(2,:)-1) - Width)/min(Height, Width);
targetPC_2D(2,:) = (Height - (2*points2d(1,:)-1))/min(Height, Width);
targetPC_2D(3,:) = ones(1, numCorres);

Shape_init.FFD = i2s_ffd_init_sym(Shape_init, Para.ffd_Res);

% Call the image-shape pair-wise aligner to obtain optimized 
[Camera_opt, Shape_FFD.ctrlPos_cur] = i2s_camera_shape_opt(targetPC_2D,...
    points3d, corresWeights, Shape_init.FFD, Camera_init, Para);
    
Shape_opt = Shape_init;

ffd_coeff_vertexPoss = i2s_ffd_basis_coeff(Shape_init.FFD,...
    double(Shape_init.vertexPoss));
Shape_opt.vertexPoss = single(Shape_init.FFD.ctrlPos_cur*...
    ffd_coeff_vertexPoss');
