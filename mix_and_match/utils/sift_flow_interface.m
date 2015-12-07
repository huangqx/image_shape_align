function [Match] = sift_flow_interface(image1, image2,...
    using_patch,... % flag indicates whether we use patch-wise corres
    patch_size,...  % the size of the patch
    patch_stride,...% stride of the patch
    max_offset) % maximum deviation of the pixel-correspondences
% 
% This function provides an interface for establishing dense pixel-wise
% correspondences from image1 to image2

if using_patch ~= 1
    % This is a simpler version --- we compute dense correspondences
    % directly on the images
    im1 = im2double(image1);
    im2 = im2double(image2);

    %figure;imshow(im1);figure;imshow(im2);
    cellsize = 3;
    gridspacing = 1;

    sift1 = mexDenseSIFT(im1, cellsize, gridspacing);
    sift2 = mexDenseSIFT(im2, cellsize, gridspacing);

    SIFTflowpara.alpha = 4*255;
    SIFTflowpara.d = 40*255;
    SIFTflowpara.gamma = 0.005*255;
    SIFTflowpara.nlevels = 4;
    SIFTflowpara.wsize = 2;
    SIFTflowpara.topwsize = 10;
    SIFTflowpara.nTopIterations = 60;
    SIFTflowpara.nIterations = 30;


    [Match.flow_12_col, Match.flow_12_row, energylist] =...
        SIFTflowc2f(sift1, sift2, SIFTflowpara);
    
    Match.flow_12_col(find(Match.flow_12_col > max_offset)) = max_offset;
    Match.flow_12_col(find(Match.flow_12_col < -max_offset)) = -max_offset;
    Match.flow_12_row(find(Match.flow_12_row > max_offset)) = max_offset;
    Match.flow_12_row(find(Match.flow_12_row < -max_offset)) = -max_offset;
else
    % perform patch sampling
    [patches, corres12] = patch_matching(...
        image1,...
        image2,...
        patch_size,...
        patch_stride,...
        max_offset);
    Match = patch_pixel_matching(image1, image2,...
        patches, corres12, patch_size);
end

function [Match] = patch_pixel_matching(image1, image2,...
    patches, corres12, patch_size)

im1 = im2double(image1);
im2 = im2double(image2);

%figure;imshow(im1);figure;imshow(im2);
cellsize = 3;
gridspacing = 1;

sift1 = mexDenseSIFT(im1, cellsize, gridspacing);
sift2 = mexDenseSIFT(im2, cellsize, gridspacing);

SIFTflowpara.alpha = 2*255;
SIFTflowpara.d = 40*255;
SIFTflowpara.gamma = 0.005*255;
SIFTflowpara.nlevels = 4;
SIFTflowpara.wsize = 2;
SIFTflowpara.topwsize = 10;
SIFTflowpara.nTopIterations = 60;
SIFTflowpara.nIterations = 30;



function [patches, corres12] = patch_matching(image1, image2,...
    patch_size, patch_stride, max_offset)

patches = mm_patch_sampling(image1, patch_size, patch_stride);
patchDess10 = mm_patch_dess(image1, patches, patch_size, 0);
patchDess11 = mm_patch_dess(image1, patches, patch_size, patch_size/2);
patchDess1 = [patchDess10; patchDess11];

patchDess20 = mm_patch_dess(image2, patches, patch_size, 0);
patchDess21 = mm_patch_dess(image2, patches, patch_size, patch_size/2);
patchDess2 = [patchDess20; patchDess21];

proj_mat = compute_proj_mat([patchDess1, patchDess2]);
patchDess1 = proj_mat*patchDess1;
patchDess2 = proj_mat*patchDess2;
patchDis = nearest_neighbor(patchDess1, patchDess2);
%
corres12 = zeros(2, size(patches, 1));
for sPId = 1:size(patches, 1)
    footId = 0;
    dis = 1e10;
    for tPId = 1:size(patches, 1)
        dcoord = patches(sPId,:) - patches(tPId,:);
        if max(dcoord) < max_offset
            if patchDis(tPId, sPId) < dis
                footId = tPId;
                dis = patchDis(tPId, sPId);
            end
        end
    end
    corres12(1, sPId) = sPId;
    corres12(2, sPId) = footId;
end

function [patchDis] = nearest_neighbor(patchDess1, patchDess2)
%
numPatches = size(patchDess1, 2);

sIds = kron(1:numPatches, ones(1, numPatches));
tIds = kron(ones(1,numPatches), 1:numPatches);
dif = patchDess1(:, sIds) - patchDess2(:, tIds);
dif = sqrt(sum(dif.*dif));
patchDis = reshape(dif, [numPatches, numPatches]);

function [proj_mat] = compute_proj_mat(patchDess)
mean_d = mean(patchDess')';
patchDess = patchDess - mean_d*single(ones(1, size(patchDess,2)));
[u,v] = eig(patchDess*patchDess');
v = diag(v);
ids = find(v > max(v)/400);
proj_mat = u(:, ids)';