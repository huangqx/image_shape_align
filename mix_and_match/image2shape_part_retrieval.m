function [Shape_init, partCorres] = image2shape_part_retrieval(...
    Image,...
    Shapes,...
    partCorresStructs,...
    Para)
% This function retrieves a subset of parts of the input shapes to
% reconstruct the input image object. This is based on solving a
% generalized set-cover problem (which ensures compatibility)
% The is a simplified (for speed concern) version of the method described
% in the Siggraph paper
% Input arguments:
%      Image: the input image
%      Shapes: the input shapes
%      partCorresStructs: the correspondence structure between each image
%                         and each shape
%      Para: please refer to the software document
% Output arguments:
%      Shape_init: the initial shape by assembling parts
%      partCorres: the dense correspondences between image pixels and the
%                  mesh points

[M_adj, IDX, setScores, Mask] = part_covering_struct(Image,...
    Shapes, partCorresStructs);

setIds = part_extraction(M_adj, IDX, setScores, Mask, partCorresStructs,...
    Para.minCoverRatio, Para.maxOverlapRatio, Para.show_imageRecons);

% Pre-align the parts
Shapes = preprocess_shape_collection(Shapes);

% Combine the retrived parts
[Shape_init, partCorres] = write_output(IDX, Shapes,...
    partCorresStructs, setIds);


function [Shape_init, partCorres] = write_output(IDX, Shapes,...
    partCorresStructs, setIds)

IDX_sub = IDX(setIds, :);
Shape_init.sym_pairs = [];
Shape_init.nonsym_ids = [];
Shape_init.faceVIds = [];
Shape_init.vertexPoss = [];
partCorres = [];
off = 0;
off_v = 0;
for i = 1:size(IDX_sub, 1)
    if IDX_sub(i,2) > 0
        shapeId = IDX_sub(i, 1);
        part1Id = IDX_sub(i, 3);
        part2Id = IDX_sub(i, 4);
        mesh1 = extract_part(Shapes{shapeId}, part1Id);
        Shape_init.meshes{off+1}.mat = mesh1.mat;
        Shape_init.meshes{off+1}.faceVIds = mesh1.faceVIds;
        Shape_init.meshes{off+1}.vertexIds = (off_v+1):(off_v+length(mesh1.vertexIds));
        partCorres{off+1}.pixels = partCorresStructs{shapeId}.partCorres{part1Id}.pixels;
        partCorres{off+1}.meshPoints = partCorresStructs{shapeId}.partCorres{part1Id}.meshPoints;
        Shape_init.faceVIds = [Shape_init.faceVIds, mesh1.faceVIds + off_v];
        Shape_init.vertexPoss = [Shape_init.vertexPoss, mesh1.vertexPoss];
        off_v = off_v + length(mesh1.vertexIds);
        
        mesh2 = extract_part(Shapes{shapeId}, part2Id);
        Shape_init.meshes{off+2}.mat = mesh2.mat;
        Shape_init.meshes{off+2}.faceVIds = mesh2.faceVIds;
        Shape_init.meshes{off+2}.vertexIds = (off_v+1):(off_v+length(mesh2.vertexIds));
        partCorres{off+2}.pixels = partCorresStructs{shapeId}.partCorres{part2Id}.pixels;
        partCorres{off+2}.meshPoints = partCorresStructs{shapeId}.partCorres{part2Id}.meshPoints;
      
        Shape_init.faceVIds = [Shape_init.faceVIds, mesh2.faceVIds + off_v];
        Shape_init.vertexPoss = [Shape_init.vertexPoss, mesh2.vertexPoss];
        Shape_init.sym_pairs = [Shape_init.sym_pairs,[off+1;off+2]];
        off_v = off_v + length(mesh2.vertexIds);
         
        off = off + 2;
    else
        shapeId = IDX_sub(i, 1);
        partId = IDX_sub(i, 3);
        mesh = extract_part(Shapes{shapeId}, partId);
        Shape_init.meshes{off+1}.mat = mesh.mat;
        Shape_init.meshes{off+1}.faceVIds = mesh.faceVIds;
        Shape_init.meshes{off+1}.vertexIds = (off_v+1):(off_v+length(mesh.vertexIds));
        partCorres{off+1}.pixels = partCorresStructs{shapeId}.partCorres{partId}.pixels;
        partCorres{off+1}.meshPoints = partCorresStructs{shapeId}.partCorres{partId}.meshPoints;
        Shape_init.faceVIds = [Shape_init.faceVIds, mesh.faceVIds + off_v];
        Shape_init.vertexPoss = [Shape_init.vertexPoss, mesh.vertexPoss];
        Shape_init.nonsym_ids = [Shape_init.nonsym_ids, off+1];
        off_v = off_v + length(mesh.vertexIds);
        off = off + 1;
    end
end
Shape_init.has_material = 1;

% Extract the parts in a greedy fashion
function [setIds] = part_extraction(M_adj,...
    IDX,...
    setScores,...
    Mask,...
    partCorresStructs,...
    minCoverRatio,...
    maxOverlapRatio,...
    show_imageRecons)
%
[Height, Width] = size(Mask);
if show_imageRecons == 1
    img_recons_r = ones(Height, Width);
    img_recons_g = ones(Height, Width);
    img_recons_b = ones(Height, Width);
end
[scores, order] = sort(setScores);
img_mask = zeros(Height, Width);
setIds = [];
rgbd = colormap('jet');

for i = 1:length(order)
    id = order(i);
    c = rgbd(max(1, min(64, floor(64*rand(1,1)))),:);
    if IDX(id, 2) > 0
        shapeId = IDX(id, 1);
        part1Id = IDX(id, 3);
        part2Id = IDX(id, 4);
        tp1 = partCorresStructs{shapeId}.partCorres{part1Id};
        tp2 = partCorresStructs{shapeId}.partCorres{part2Id};
        if length(tp2.pixels2) == 0 && length(tp2.pixels2) == 0
            continue;
        end
        pixels = [tp1.pixels, tp2.pixels];
        pixelsId = (pixels(2,:)-1)*Height + pixels(1,:);
        pixels2 = [tp1.pixels2, tp2.pixels2];
        pixels2Id = (pixels2(2,:)-1)*Height + pixels2(1,:);
        ratio = length(pixelsId)/length(pixels2Id);
        if ratio < minCoverRatio
            continue;
        end
        ratio2 = sum(img_mask(pixels2Id))/length(pixels2Id);
        if ratio2 > maxOverlapRatio
            continue;
        end
        if show_imageRecons == 1
            img_recons_r(pixels2Id) = c(1);
            img_recons_g(pixels2Id) = c(2);
            img_recons_b(pixels2Id) = c(3);
        end
        img_mask(pixels2Id) = 1;
        setIds = [setIds, id];
    else
        shapeId = IDX(id, 1);
        part1Id = IDX(id, 3);
        tp1 = partCorresStructs{shapeId}.partCorres{part1Id};
        if length(tp1.pixels2) == 0
            continue;
        end
        pixels = tp1.pixels;
        pixelsId = (pixels(2,:)-1)*Height + pixels(1,:);
        pixels2 = tp1.pixels2;
        pixels2Id = (pixels2(2,:)-1)*Height + pixels2(1,:);
        ratio = length(pixelsId)/length(pixels2Id);
        if ratio < minCoverRatio
            continue;
        end
        ratio2 = sum(img_mask(pixels2Id))/length(pixels2Id);
        if ratio2 > maxOverlapRatio
            continue;
        end
        if show_imageRecons == 1
            img_recons_r(pixels2Id) = c(1);
            img_recons_g(pixels2Id) = c(2);
            img_recons_b(pixels2Id) = c(3);
        end
        img_mask(pixels2Id) = 1;
        setIds = [setIds, id];
    end
end
img_recons = zeros(Height, Width, 3);
img_recons(:,:,1) = img_recons_r;
img_recons(:,:,2) = img_recons_g;
img_recons(:,:,3) = img_recons_b;
if show_imageRecons == 1
    figure(2);
    imshow(img_recons);
end

function [M_adj, IDX, setScores, Mask] = part_covering_struct(Image,...
    Shapes, partCorresStructs)
% Determine the sets, which consists of non-symmetric parts and symmetric
% part pairs
% Input arguments:
%     Image: the input image
%     Shapes: the input shapes
%     partCorresStructs: the dense correspondences
% Output arguments:
%     M_adj: the adjacency information between the sets of the pixels
%     IDX: the index structure between the sets of the parts
%     setScores: the scores of the sets
%     Mask: the mask associated with the input image object

% The number of sets is equivalent to the number of symmetric part pairs
% and the number of non-symmetric parts
[Height, Width, k] = size(Image.im);
numSets = 0;
for i = 1 : length(partCorresStructs)
    numSets = numSets + size(Shapes{i}.sym_pairs, 2);
    numSets = numSets + length(Shapes{i}.nonsym_ids);
end


% Compute the propagated mask for the input image object
Mask = zeros(Height, Width);
for id = 1:length(partCorresStructs)
    pcs = partCorresStructs{id};
    for j = 1:length(pcs.partCorres)
        pc = pcs.partCorres{j};
        if length(pc.pixels) > 0
            pixelIds = (pc.pixels(2,:)-1)*Height + pc.pixels(1,:);
            Mask(pixelIds) = 1;
        end
    end
end

numPixels = Height*Width;

% The association matrix between the sets and the image object
M_adj = sparse(numSets, numPixels);

setScores = zeros(1, numSets); %
IDX = zeros(numSets, 4); % a index table for the sets

setId = 0;
for shapeId = 1 : length(partCorresStructs)
    % Go through symmetric parts and non symmetric parts
    for j = 1 : size(Shapes{shapeId}.sym_pairs, 2)
        setId = setId + 1;
        part1Id = Shapes{shapeId}.sym_pairs(1, j);
        part2Id = Shapes{shapeId}.sym_pairs(2, j);
        tp1 = partCorresStructs{shapeId}.partCorres{part1Id};
        tp2 = partCorresStructs{shapeId}.partCorres{part2Id};
        
        % Set the index table
        IDX(setId, 1) = shapeId;
        IDX(setId, 2) = 1;
        IDX(setId, 3) = Shapes{shapeId}.sym_pairs(1, j);
        IDX(setId, 4) = Shapes{shapeId}.sym_pairs(2, j);
        pixels = [tp1.pixels, tp2.pixels];
        %
        if length(pixels) > 0
            pixels2 = [tp1.pixels2, tp2.pixels2];
            pixel2Ids = (pixels2(2,:)-1)*Height + pixels2(1,:);
            M_adj(setId, pixel2Ids) = 1;
            % Compute the set score
            % Compute the similarity score
            simDist = 1;
            if isfield(tp1, 'similarityDist') == 1 &&...
                    isfield(tp2, 'similarityDist') == 0
                simDist = tp1.similarityDist;
            end
            if isfield(tp1, 'similarityDist') == 0 &&...
                    isfield(tp2, 'similarityDist') == 1
                simDist = tp2.similarityDist;
            end
            if isfield(tp1, 'similarityDist') == 1.&&...
                    isfield(tp2, 'similarityDist') == 1
                simDist = max(tp1.similarityDist, tp2.similarityDist);
            end
            % the overlapping ratio
            overlapRatio = size(pixels,2)/size(pixels2, 2);
            score = simDist*power(1/overlapRatio, 2)/size(pixels,2);
            setScores(setId) = score;
        end
    end
    for j = 1 : length(Shapes{shapeId}.nonsym_ids)
        setId = setId + 1;
        partId = Shapes{shapeId}.nonsym_ids(j);
        tp0 = partCorresStructs{shapeId}.partCorres{partId};
        
        % Set the index table
        IDX(setId, 1) = shapeId;
        IDX(setId, 2) = -1;
        IDX(setId, 3) = Shapes{shapeId}.nonsym_ids(j);
        IDX(setId, 4) = 0;
        %
        pixels = tp0.pixels;
        if length(pixels) > 0
            pixels2 = tp0.pixels2;
            pixel2Ids = (pixels2(2,:)-1)*Height + pixels2(1,:);
            M_adj(setId, pixel2Ids) = 1;
            % Compute the set score
            % Compute the similarity score
            % the overlapping ratio
            overlapRatio = size(pixels,2)/size(pixels2, 2);
            score = tp0.similarityDist*power(1/overlapRatio, 2)/size(pixels,2);
            setScores(setId) = score;
        end
    end
end

ids = find(setScores == 0);
setScores(ids) = 1e3;

function [mesh] = extract_part(Shape, partId)

mesh = Shape.meshes{partId};
mesh.vertexPoss = Shape.vertexPoss(:, mesh.vertexIds);