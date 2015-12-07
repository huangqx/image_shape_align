function [patchDess] = mm_patch_dess(img, patches, patch_size, pad)
% expand the image
[m, n, k] = size(img);
m_exp = m + 2*pad;
n_exp = n + 2*pad;

patch_size2 = patch_size + 2*pad;

img_exp = 255 - uint8(zeros(m_exp, n_exp, 3));
img_exp((pad+1):(pad+m),(pad+1):(pad+n),:) = img;

patchDess = single(zeros(1024, size(patches,1)));

for pId = 1 : size(patches, 1)
    rowIds = patches(pId,1):(patches(pId,1) + patch_size2 - 1);
    colIds = patches(pId,2):(patches(pId,2) + patch_size2 - 1);

    patch = img_exp(rowIds, colIds, :);
    patch = imResample(single(patch), [128, 128])/255;
    H = hog(patch, 32, 16);    
    H = reshape(H, [1024,1]);
    patchDess(:, pId) = H; 
end