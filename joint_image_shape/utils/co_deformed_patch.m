function [img] = co_deformed_patch(patch, ffd)

patch = im2double(patch);
[m, n, k] = size(patch);

pts = ffd.ctrlPts*ffd.basis';

patch_r = reshape(patch(:,:,1), [1, size(pts,2)]);
patch_g = reshape(patch(:,:,2), [1, size(pts,2)]);
patch_b = reshape(patch(:,:,3), [1, size(pts,2)]);

lattice = [pts; patch_r;patch_g;patch_b];

para.lowerCorner = [-8,-8];
[para.latticeDimY, para.latticeDimX, k] = size(patch);
para.nCols = para.latticeDimX + 16;
para.nRows = para.latticeDimY + 16;
para.pixelGap = 1;
img = lattice_2_image(lattice, para);

