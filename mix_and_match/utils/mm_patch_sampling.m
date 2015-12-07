function [patches] = mm_patch_sampling(img, patch_size, patch_stride)
%
[height, width, k] = size(img);

nRows = floor((height - patch_size)/patch_stride)+1;
nCols = floor((width - patch_size)/patch_stride)+1;

rowOffs = [0:patch_stride:(patch_stride*(nRows-1))];
colOffs = [0:patch_stride:(patch_stride*(nCols-1))];

if abs(rowOffs(nRows) - height+patch_size) < patch_stride/2
    rowOffs(nRows) = height - patch_size;
else
    rowOffs = [rowOffs, height - patch_size];
end

if abs(colOffs(nCols) - width+patch_size) < patch_stride/2
    colOffs(nCols) = width - patch_size;
else
    colOffs = [colOffs, width - patch_size];
end

dim = length(rowOffs)*length(colOffs);

off = 0;
patches = zeros(dim, 2);
for i = 1:length(rowOffs)
    for j = 1:length(colOffs)
        off = off + 1;
        patches(off, :) = [rowOffs(i), colOffs(j)];
    end
end
patches = patches + 1;