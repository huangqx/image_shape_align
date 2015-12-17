function [Clusters] = jis_view_clustering(ImageCameras, viewDirectRadius)
% Divide input images into clusters with similar view positions.
% Input arguments:
%      ImageCameras: the camera configurations associated with input images
%      directDis: the threshold on

numCameras = length(ImageCameras);
viewDirections = zeros(3, numCameras);

for camId = 1:numCameras
    vec = ImageCameras{camId}.origin - ImageCameras{camId}.lookAt;
    vec = vec/norm(vec);
    viewDirections(:, camId) = vec;
end

d = viewDirections - viewDirections(:,1)*ones(1, numCameras);
[s, footId] = max(sum(d.*d));

IDX = zeros(1, numCameras);
dis = ones(1, numCameras)*1e10;
footIds = [footId];
for i = 1 : numCameras
    d = viewDirections - viewDirections(:,footId)*ones(1, numCameras);
    d = sum(d.*d);
    ids = find(d < dis);
    IDX(ids) = i;
    dis(ids) = d(ids);
    [s, footId] = max(dis);
    footIds = [footIds, footId];
    if s < viewDirectRadius
        break;
    end
end

numClusters = max(IDX);
for id = 1 : numClusters
    ids = find(IDX == id);
    Clusters{id}.imageIds = ids;
    Clusters{id}.Camera = ImageCameras{footIds(id)};
    % adjust the camera so that the rendered image enclose the image object
    Clusters{id}.Camera.scale = Clusters{id}.Camera.scale*1.1;
end


