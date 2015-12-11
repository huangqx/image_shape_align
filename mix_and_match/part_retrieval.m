function [Shape_init, partCorres] = part_retrieval(Image, Shapes,...
    partCorresStructs, Para)
% This function retrieves a subset of parts of the input shapes to
% reconstruct the input image object. This is based on solving a
% generalized set-cover problem

[M_adj, IDX, setScores, Mask] = part_covering_struct(Image,...
    Shapes, partCorresStructs);

setIds = part_extraction(M_adj, IDX, setScores, Mask, partCorresStructs);

h = 10;

% Extract the parts in a greedy fashion
function [setIds] = part_extraction(M_adj,...
    IDX,...
    setScores,...
    Mask,...
    partCorresStructs)
%
[Height, Width] = size(Mask);
img_recons = ones(Height, Width, 3);
img_mask = zeros(Height, Width);
[scores, order] = sort(setScores);

setIds = [];
for i = 1:length(order)
    id = order(i);
    if IDX(id, 2) > 0
        shapeId = IDX(id, 1);
        part1Id = IDX(id, 3);
        part2Id = IDX(id, 4);
        tp1 = partCorresStructs{shapeId}.partCorres{part1Id};
        tp2 = partCorresStructs{shapeId}.partCorres{part2Id};
        if isfield(tp1, 'pixels2') == 0 || isfield(tp2, 'pixels2') == 0
            continue;
        end
        pixels = [tp1.pixels, tp2.pixels];
        pixelsId = (pixels(2,:)-1)*Height + pixels(1,:);
        pixels2 = [tp1.pixels2, tp2.pixels2];
        pixels2Id = (pixels2(2,:)-1)*Height + pixels2(1,:);
        ratio = length(pixelsId)/length(pixels2Id);
        if ratio < 0.75
            continue;
        end
        ratio2 = sum(img_mask(pixels2Id))/length(pixels2Id);
        if ratio2 > 0.15
            continue;
        end
        img_mask(pixels2Id) = 1;
        setIds = [setIds, id];
    else
        shapeId = IDX(id, 1);
        part1Id = IDX(id, 3);
        tp1 = partCorresStructs{shapeId}.partCorres{part1Id};
        if isfield(tp1, 'pixels2') == 0
            continue;
        end
        pixels = tp1.pixels;
        pixelsId = (pixels(2,:)-1)*Height + pixels(1,:);
        pixels2 = tp1.pixels2;
        pixels2Id = (pixels2(2,:)-1)*Height + pixels2(1,:);
        ratio = length(pixelsId)/length(pixels2Id);
        if ratio < 0.75
            continue;
        end
        ratio2 = sum(img_mask(pixels2Id))/length(pixels2Id);
        if ratio2 > 0.15
            continue;
        end
        img_mask(pixels2Id) = 1;
        setIds = [setIds, id];
    end
end
h = 10;

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
        IDX(setId, 1) = i;
        IDX(setId, 2) = -1;
        IDX(setId, 3) = Shapes{i}.nonsym_ids(j);
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