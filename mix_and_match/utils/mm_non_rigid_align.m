function [Shape_opt, Camera_opt] = mm_non_rigid_align(Image, Shape_init,...
    Camera_init, partCorres)
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
pt2d = zeros(2, numCorres);
pt3d = zeros(3, numCorres);

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
    pt2d(:, (numCorres+1):(numCorres+nf)) = pc.pixels;
    pt3d(:, (numCorres+1):(numCorres+nf)) = p1.*t1 + p2.*t2 + p3.*t3;
    numCorres = numCorres + nf;
end

h = 10;
