function [partCorresStruct] = image2shape_dense_corres(Image,... 
    Shape, Camera, Para)
% Compute part-aware dense correspondences between the 3D shapes and the 
% input image object
% Input arguments:
%   Image:   The input image
%   Shape:   The input shapes described in a world coordinate system
%   Camera:  The corresponding camera configurations obtained either from
%            image-shape alignment or pose estimation
%   Para:    Please refer to the software document
% Output argment:
%   partCorresStruct:  Compute part-wise correspondences between each shape
%                      part and the corresponding pixels in images


if 1
    % Color different parts with different colors
    Shape.has_material = 1;
end

[renderImage, pixelPartIds2, meshPoints] = mm_render_shape(Shape, Camera);

% Compute patch-wise sift-flow
Match = sift_flow_interface(Image.im, renderImage, 0,...
    Para.patch_size, Para.patch_stride, Para.max_offset);


[Height, Width, k] = size(renderImage);
ROWId = kron(ones(1, Width), 1:Height);
COLId = kron(1:Width, ones(1, Height));
dRow = reshape(Match.flow_12_col, [1, Height*Width]);
dCol = reshape(Match.flow_12_row, [1, Height*Width]);
ROWId = ROWId + dRow;
COLId = COLId + dCol;
ids = find(ROWId <= Height & ROWId >= 1 & COLId <= Width & COLId >= 1);

pixelIds = (COLId(ids)-1)*Height + ROWId(ids);
pixelPartIds = zeros(Height, Width);
pixelPartIds(ids) = pixelPartIds2(pixelIds);

% apply the mask on the source image
mask = double(rgb2gray(Image.im) < 240);
pixelPartIds = pixelPartIds.*mask;

h = 10;

if 1 % used for debugging
    %
    if 1
        hold on;
        subplot(1,3,1);
        imshow(Image.im);
        subplot(1,3,2);
        imagesc(pixelPartIds);
        daspect([1,1,1]);
        subplot(1,3,3);
        imagesc(pixelPartIds2);
        daspect([1,1,1]);
    end
end
partCorresStruct = [];
%
%

%




%figure;imshow(im1);figure;imshow(im2);


% Compute sift-flow correspondences
