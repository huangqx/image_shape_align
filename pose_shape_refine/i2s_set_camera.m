function [Camera_out] = i2s_set_camera(image, Camera_in, rHeights, rWidths)
% Set the camera window size in Matlab. Note that we obtain an image
% of the same size, one has to draw a bigger image

[nHeight_inner, nWidth_inner, k] = size(image);

Camera_out = Camera_in;
Camera_out.nHeight_inner = nHeight_inner;
Camera_out.nWidth_inner = nWidth_inner;

tp = find(rHeights == nHeight_inner);
Camera_out.nHeight = tp(1);

tp = find(rWidths == nWidth_inner);
Camera_out.nWidth = tp(1);
