function [ffds, hogDess] = dense_multi_image_corres(images, Para)
% Align the input images in a common coordinate system so that the
% corresponding features are aligned
% Input argument:
%    images: the set of images consists of both natural and rendered images
% Output argument:
%    ffds: the free-form deformations that align all the input images

% Compute pair-wise dense correspondences between pairs of images
fprintf('  pairwise matching...\n');
[maps, hogDess] = co_pairwise_matching(images,...
    Para.knn_multi_image_matching);
fprintf('  joint matching...\n');
ffds = co_ffd_joint_align(images, maps, Para);