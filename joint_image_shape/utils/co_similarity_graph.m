function [G] = co_similarity_graph(dess, knn)
% This function computes a similarity graph among the input images by
% utilizing the image descriptors
knn = min(size(dess,2)-1, knn);
numImages = size(dess, 2);
sIds = kron(1:numImages, ones(1, numImages));
tIds = kron(ones(1, numImages), 1:numImages);
simMat = dess(:, sIds) - dess(:, tIds);
simMat = sum(simMat.*simMat);
simMat = reshape(simMat, [numImages, numImages]);

G = sparse(numImages, numImages);
for sId = 1 : numImages
    [s, tIds] = sort(simMat(sId, :));
    tIds = tIds(2:(knn+1));
    G(sId, tIds) = s(2:(knn+1));
end
G = max(G, G');
[rows, cols, vals] = find(G);
vals = exp(-vals/2/median(vals));
G = sparse(rows, cols, vals);

