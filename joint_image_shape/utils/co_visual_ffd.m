function [imgs] = co_visual_ffd(patches, ffds)

% View patch
figure(1)
subplot(2,3,1)
imshow(patches{1});
subplot(2,3,2)
imshow(patches{2});
subplot(2,3,3)
imshow(patches{3});
subplot(2,3,4)
imshow(patches{4});
subplot(2,3,5)
imshow(patches{5});
subplot(2,3,6)
imshow(patches{6});

% subplot(2,4,7)
% imshow(patches{7});
% subplot(2,4,8)
% imshow(patches{8});

figure(2)
subplot(2,3,1)
img = deformed_patch(patches{1}, ffds{1});
imshow(img);
imgs{1} = img;

subplot(2,3,2)
img = deformed_patch(patches{2}, ffds{2});
imshow(img);
imgs{2} = img;

subplot(2,3,3)
img = deformed_patch(patches{3}, ffds{3});
imshow(img);
imgs{3} = img;

subplot(2,3,4)
img = deformed_patch(patches{4}, ffds{4});
imshow(img);
imgs{4} = img;

subplot(2,3,5)
img = deformed_patch(patches{5}, ffds{5});
imshow(img);
imgs{5} = img;

subplot(2,3,6)
img = deformed_patch(patches{6}, ffds{6});
imshow(img);
imgs{6} = img;

%tp = (double(imgs{1}) + double(imgs{2}) + double(imgs{3})...
%    + double(imgs{4}) + double(imgs{5})+ double(imgs{6}))/6;
%imgs{7} = tp/max(max(max(tp)));

% img = deformed_patch(patches{5}, ffds{5});
% imshow(img);
% imgs{6} = img;
% 
% subplot(2,4,7)
% img = deformed_patch(patches{7}, ffds{7});
% imshow(img);
% imgs{7} = img;
% 
% subplot(2,4,8)
% img = deformed_patch(patches{8}, ffds{8});
% imshow(img);
% imgs{8} = img;


function [img] = deformed_patch(patch, ffd)

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

% map = maps{id};
% 
% patch_s = patches{map.sId};
% patch_t = patches{map.tId};
% 
% patch_sft = 255-uint8(zeros(size(patch_s)));
% [m,n,k] = size(patch_sft);
% 
% for col = 1 : n
%     for row = 1 : m
%         col1 = map.ts_col(row, col) + col;
%         row1 = map.ts_row(row, col) + row;
%         if 1 <= col1 && col1 <= n && 1 <= row1 && row1 <= m
%             patch_sft(row1, col1,:) = patch_t(row, col,:);
%         end
%     end
% end
% 
% patch_tfs = 255-uint8(zeros(size(patch_t)));
% [m,n,k] = size(patch_tfs);
% 
% for col = 1 : n
%     for row = 1 : m
%         col1 = map.st_col(row, col) + col;
%         row1 = map.st_row(row, col) + row;
%         if 1 <= col1 && col1 <= n && 1 <= row1 && row1 <= m
%             patch_tfs(row1, col1,:) = patch_s(row, col,:);
%         end
%     end
% end
% 
% figure(1);
% subplot(2,2,1);
% imshow(patch_s);
% 
% subplot(2,2,2);
% imshow(patch_sft);
% 
% subplot(2,2,3);
% imshow(patch_t);
% 
% subplot(2,2,4);
% imshow(patch_tfs);