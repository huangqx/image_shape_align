function [Shape_init, partCorres] = part_retrieval(Image, Shapes,...
    partCorresStructs, Para)
% This function retrieves a subset of parts of the input shapes to
% reconstruct the input image object. This is based on solving a
% generalized set-cover problem

numSets = 0;
for i = 1 : length(partCorresStructs)
    numSets = numSets + size(Shapes{i}.sym_pairs, 2);
    numSets = numSets + length(Shapes{i}.nonsym_ids);
end

[Height, Width, k] = size(Image.im);

numPixels = Height*Width;
Adj = sparse(numSets, numPixels);

scores = ones(1, numSets)*10;

off = 0;
for i = 1 : length(partCorresStructs)
    for j = 1 : size(Shapes{i}.sym_pairs, 2)
        off = off + 1;
        tp1 = partCorresStructs{i}.partCorres{Shapes{i}.sym_pairs(1, j)};
        tp2 = partCorresStructs{i}.partCorres{Shapes{i}.sym_pairs(2, j)};
        pixels2 = [tp1.pixels, tp2.pixels];
        if length(pixels2) > 0
            pixelIds = (pixels2(2,:)-1)*Height + pixels2(1,:);
            Adj(off, pixelIds) = 1;
            if isfield(tp1, 'similarityDist') == 1
                scores(off) = min(scores(off), tp1.similarityDist);
            end
            if isfield(tp2, 'similarityDist') == 1
                scores(off) = min(scores(off), tp2.similarityDist);
            end
        end
    end
    for j = 1 : length(Shapes{i}.nonsym_ids)
        off = off + 1;
        tp0 = partCorresStructs{i}.partCorres{Shapes{i}.nonsym_ids(j)};
        pixels2 = tp0.pixels;
        if length(pixels2) > 0
            pixelIds = (pixels2(2,:)-1)*Height + pixels2(1,:);
            Adj(off, pixelIds) = 1;
            if isfield(tp0, 'similarityDist') == 1
                scores(off) = min(scores(off), tp0.similarityDist);
            end
        end
    end
end

h = 10;