function [img] = lattice_2_image(lattice, para)

paras = [para.latticeDimX, para.latticeDimY, para.lowerCorner(1), para.lowerCorner(2), para.nCols, para.nRows, para.pixelGap];
img1 = lattice2image(lattice, paras);

img = zeros(para.nRows, para.nCols, 3);
img(:,:,1) = reshape(img1(1,:), [para.nRows, para.nCols]);
img(:,:,2) = reshape(img1(2,:), [para.nRows, para.nCols]);
img(:,:,3) = reshape(img1(3,:), [para.nRows, para.nCols]);
