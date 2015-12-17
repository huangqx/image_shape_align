function [corres] = establish_2d_3d_corres(ffd1, ffd2, scan2)
% Estimate dense correspondences between im1 and rgbd scan
% Input arguments:
%  ffd1: the free-from deformation of image 1
%  ffd2: the free-from deformation of image 2
%  scan2: the point cloud that is in correspondence with im2
% output argument:
%  corres.pixelIds: the image pixel indices (column-major)
%  corres: meshPoints: the corresponding mesh points ([faceIds; t1; t2])

% the deformed pixel coordinates
pt1 = ffd1.ctrlPts*ffd1.basis';
pt2 = ffd2.ctrlPts*ffd2.basis';

% find valid pixels of the rendered images
ids = find(scan2(4,:) < 1e5);
pt2 = pt2(:, ids);
meshPoints = scan2(1:3, ids);

% perform nearest-neighbor search
[IDX_12, D_12] = knnsearch(pt2', pt1');
corres.pixelIds = find(D_12 < 2)';
corres.meshPoints = meshPoints(:, IDX_12(corres.pixelIds));

