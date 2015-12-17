function [maps, dess] = co_pairwise_matching(patches, knn)

% Parameters
SIFTflowpara.alpha=2*255;
SIFTflowpara.d=20*255;
SIFTflowpara.gamma=0.02*255;
SIFTflowpara.nlevels=4;
SIFTflowpara.wsize=2;
SIFTflowpara.topwsize=10;
SIFTflowpara.nTopIterations = 60;
SIFTflowpara.nIterations= 30;
cellsize=3;
gridspacing=1;

% ims
numPatches = length(patches);
dess = zeros(1024, numPatches);
for i = 1:length(patches)
    dess(:, i) = co_hog_des(patches{i});
end
G = co_similarity_graph(dess, knn);
[rows, cols, vals] = find(G);
edges = [rows, cols, vals]';
edges = edges(:, find(rows < cols));

[scales, scaled_patches] = normalize_patches(patches);

sifts = cell(1, length(scaled_patches));
for i = 1 : length(scaled_patches)
    patch = scaled_patches{i};
    sifts{i} = mexDenseSIFT(im2double(patch), cellsize,gridspacing);
end

maps = cell(1, size(edges, 2));
for id = 1 : size(edges, 2)
    maps{id}.sId = edges(1, id);
    maps{id}.tId = edges(2, id);
    maps{id}.weight = edges(3, id);
    maps{id}.corres_st = sift_corres(sifts{maps{id}.sId},...
        sifts{maps{id}.tId},...
        scales(:, maps{id}.sId),...
        scales(:, maps{id}.tId),...
        SIFTflowpara);
    maps{id}.corres_ts = sift_corres(sifts{maps{id}.tId},...
        sifts{maps{id}.sId},...
        scales(:, maps{id}.tId),...
        scales(:, maps{id}.sId),...
        SIFTflowpara);
%    fprintf('[%d, %d]\n', id, size(edges, 2));
end

function [corres] = sift_corres(sift_s, sift_t, scale_s, scale_t, SIFTflowpara)

[dif_col, dif_row, energylist] = SIFTflowc2f(sift_s, sift_t, SIFTflowpara);
    
dif_col(find(dif_col > 128)) = 128;
dif_col(find(dif_col < -128)) = -128;
dif_row(find(dif_row > 128)) = 128;
dif_row(find(dif_row < -128)) = -128;
    
[ms, ns, k] = size(sift_s);
[mt, nt, k] = size(sift_t);
row_s = kron((1:ms)', ones(1,ns));
col_s = kron(ones(ms,1), 1:ns);
row_t = row_s + dif_row;
col_t = col_s + dif_col;

row_s = reshape(row_s, [1, ms*ns]);
col_s = reshape(col_s, [1, ms*ns]);
row_t = reshape(row_t, [1, ms*ns]);
col_t = reshape(col_t, [1, ms*ns]);

tIds = (col_t-1)*mt + row_t;
ids = find(tIds <= mt*nt & tIds >= 1);
ids = ids(1:8:length(ids));
row_s = row_s(ids);
row_t = row_t(ids);
col_s = col_s(ids);
col_t = col_t(ids);

row_s = (row_s - 0.5)/scale_s(1);
col_s = (col_s - 0.5)/scale_s(2);

row_t = (row_t - 0.5)/scale_t(1);
col_t = (col_t - 0.5)/scale_t(2);

row_s1 = max(1, min(ms, floor(row_s +0.5)));
col_s1 = max(1, min(ns, floor(col_s +0.5)));
row_t1 = max(1, min(mt, floor(row_t +0.5)));
col_t1 = max(1, min(nt, floor(col_t +0.5)));

ms = floor((ms/scale_s(1))+0.5);
mt = floor((mt/scale_t(1))+0.5);

sIds = (col_s1 - 1)*ms + row_s1;
tIds = (col_t1 - 1)*mt + row_t1;
corres = [sIds; tIds; row_s; col_s; row_t; col_t];

function [scales, scaled_patches] = normalize_patches(patches)
%
numPatches = length(patches);
bboxes = zeros(2, numPatches);
for i = 1 : numPatches
    [m,n,k] = size(patches{i});
    scale = 500/norm([m,n]);
    bboxes(1, i) = m*scale;
    bboxes(2, i) = n*scale;
end
%
mean_row = floor(mean(bboxes(1,:)));
mean_col = floor(mean(bboxes(2,:)));

scales = zeros(2, numPatches);
scaled_patches = cell(1, numPatches);
for i = 1 : numPatches
    [m,n,k] = size(patches{i});
    scales(1,i) = mean_row/m;
    scales(2,i) = mean_col/n;
    scaled_patches{i} = imresize(patches{i}, [mean_row, mean_col]);
end
