function [Match] = image2shape_dense_corres(Image,... 
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

[renderImage_aligned, pixelPartIds] = deform_image(renderImage,...
    pixelPartIds2, Match);

if 1
    tp = (im2double(renderImage_aligned) + im2double(Image.im))/2;
    imshow(tp);
end

nf = size(Shape.faceVIds, 2);
facePartIds = zeros(1, nf);
nf= 0;
offs_part = zeros(1, length(Shape.meshes));
for partId = 1:length(Shape.meshes)
    mesh = Shape.meshes{partId};
    nf_mesh = size(mesh.faceVIds, 2);
    facePartIds((nf+1):(nf+nf_mesh)) = partId;
    offs_part(partId) = nf;
    nf = nf + nf_mesh;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract part-wise correspondences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the mask on the source image
mask = double(rgb2gray(Image.im) < 240);
pixelPartIds0 = pixelPartIds.*mask;

numParts = size(Shape.sym_pairs,2)*2 + length(Shape.nonsym_ids);
[Height, Width] = size(mask);
for partId = 1:numParts
    [rows, cols, vals] = find(pixelPartIds0 == partId);
    Match.partCorres{partId}.pixels = [];
    Match.partCorres{partId}.meshPoints = [];
    Match.partCorres{partId}.pixels2 = [];
    boxes{partId} = [];
    if length(rows) > 0
        pixelIds = (cols'-1)*Height + rows';
        Match.partCorres{partId}.pixels = [rows, cols]';
        [rows2, cols2, vals] = find(pixelPartIds == partId);
        Match.partCorres{partId}.pixels2 = [rows2, cols2]';
        rows = rows + Match.flow_12_row(pixelIds)';
        cols = cols + Match.flow_12_col(pixelIds)';
        pixelIds2 = (cols'-1)*Height + rows';
        Match.partCorres{partId}.meshPoints = meshPoints(:, pixelIds2);
        fids_global = Match.partCorres{partId}.meshPoints(1,:);
        partId_local = facePartIds(fids_global);

        Match.partCorres{partId}.meshPoints(1,:) = fids_global...
            - offs_part(partId_local(1));
        boxes{partId} = [min(cols), min(rows), max(cols), max(rows)];
    end
end

simiDists = compare_parts(Image.im, renderImage_aligned, boxes,...
    Para.patch_pad);

for partId = 1:numParts
    if length(Match.partCorres{partId}.pixels) ~= 0
        Match.partCorres{partId}.similarityDist = simiDists(partId);
    end
end
Match.renderImage_aligned = renderImage_aligned;

if 0 % used for debugging
    %
    hold on;
    subplot(1,3,1);
    imshow(Image.im);
    subplot(1,3,2);
    imagesc(pixelPartIds0);
    daspect([1,1,1]);
    subplot(1,3,3);
    imagesc(pixelPartIds2);
    daspect([1,1,1]);
end

function [simiDists] = compare_parts(inputImage, renderedImage, boxes, pad)
%
[m,n,k] = size(inputImage);
m1 = m + 2*pad;
n1 = n + 2*pad;

image1 = 255 - uint8(zeros(m1, n1, 3));
image1((pad+1):(pad+m),(pad+1):(pad+n),:) = inputImage;

image2 = 255 - uint8(zeros(m1, n1, 3));
image2((pad+1):(pad+m),(pad+1):(pad+n),:) = renderedImage;

simiDists = zeros(1, length(boxes));
for i = 1:length(boxes)
    if length(boxes{i}) ~= 0
        box = boxes{i};
        colIds = box(1):(box(3) + 2*pad);
        rowIds = box(2):(box(4) + 2*pad);
        patch1 = image1(rowIds, colIds, :);
        patch2 = image2(rowIds, colIds, :);
        
        s = length(rowIds)/length(colIds);
        if s < 1
            nr = max(2, floor(4*s + 0.5));
            nc = 4;
        else
            nr = 4;
            nc = max(2, floor(4/s + 0.5));
        end
        patch1 = imResample(single(patch1), [64*nr, 64*nc])/255;
        H1 = hog(patch1, 64, 16);    
        H1 = reshape(H1, [64*nc*nr,1]);
        
        patch2 = imResample(single(patch2), [64*nr, 64*nc])/255;
        H2 = hog(patch2, 64, 16);    
        H2 = reshape(H2, [64*nc*nr,1]);
        simiDists(i) = norm(H1 - H2)/sqrt(nc*nr);
    end
end

function [renderImage_aligned, pixelPartIds] = deform_image(renderImage,...
    pixelPartIds2, Match)

[Height, Width, k] = size(renderImage);
ROWId = kron(ones(1, Width), 1:Height);
COLId = kron(1:Width, ones(1, Height));
dRow = reshape(Match.flow_12_row, [1, Height*Width]);
dCol = reshape(Match.flow_12_col, [1, Height*Width]);
ROWId = ROWId + dRow;
COLId = COLId + dCol;
ids = find(ROWId <= Height & ROWId >= 1 & COLId <= Width & COLId >= 1);
pixelIds = (COLId(ids)-1)*Height + ROWId(ids);
pixelPartIds = zeros(Height, Width);
pixelPartIds(ids) = pixelPartIds2(pixelIds);

renderImage_aligned_r = 255 - uint8(zeros(Height, Width));
renderImage_aligned_g = 255 - uint8(zeros(Height, Width));
renderImage_aligned_b = 255 - uint8(zeros(Height, Width));
tp = renderImage(:,:,1);
renderImage_aligned_r(ids) = tp(pixelIds);
tp = renderImage(:,:,2);
renderImage_aligned_g(ids) = tp(pixelIds);
tp = renderImage(:,:,3);
renderImage_aligned_b(ids) = tp(pixelIds);
renderImage_aligned = 255 - uint8(zeros(Height, Width, 3));
renderImage_aligned(:,:,1) = renderImage_aligned_r;
renderImage_aligned(:,:,2) = renderImage_aligned_g;
renderImage_aligned(:,:,3) = renderImage_aligned_b;
