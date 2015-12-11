function [meshes] = mm_parts_from_shape(Shape, err)
%
numVertices = size(Shape.vertexPoss, 2);
ids = zeros(1, numVertices);

id = 1;
for i = 1:numVertices
    if ids(i) > 0
        continue;
    end
    d = Shape.vertexPoss(:, i)*single(ones(1, numVertices)) - Shape.vertexPoss;
    tp = find(sqrt(sum(d.*d)) < err);
    if length(tp) >= 1
        tp = tp(find(tp >= i));
        ids(tp) = id;
    end
    id = id + 1;
end

numV = id - 1;
fVIds = ids(Shape.faceVIds);
rows = [fVIds(1,:), fVIds(2,:), fVIds(3,:)];
ools = [fVIds(2,:), fVIds(3,:), fVIds(1,:)]; 
G = sparse(rows, ools, ones(1, length(rows)), numV,numV);
G = max(G, G');
IDX = dis_components(G);

facePartIds = IDX(fVIds(1,:));
numParts = max(facePartIds);
for i = 1:numParts
    fvIds = Shape.faceVIds(:, find(facePartIds == i));
    flags = zeros(1, numVertices);
    flags(fvIds) = 1;
    mesh.vertexIds = find(flags == 1);
    flags(mesh.vertexIds) = 1:length(mesh.vertexIds);
    mesh.faceVIds = flags(fvIds);
    mesh.mat = Shape.meshes{1}.mat;
    meshes{i} = mesh;
end

function [IDX] = dis_components(G)

numV = size(G, 1);
IDX = zeros(1, numV);
flags = zeros(1, numV);
headId = 1;

clusterId = 0;
while 1
    while headId <= numV
        if flags(headId) == 0
            break;
        end
        headId = headId + 1;
    end
    if headId > numV
        break;
    end
    clusterId = clusterId + 1;
    clus = [headId];
    flags(headId) = 1;
    fringe_start = 1;
    fringe_end = 1;
    while 1
        for i = fringe_start:fringe_end
            vId = clus(i);
            nIds = find(G(vId, :));
            for j = 1:length(nIds)
                nId = nIds(j);
                if flags(nId) == 0
                    clus = [clus, nId];
                    flags(nId) = 1;
                end
            end
        end
        fringe_start = fringe_end + 1;
        fringe_end = length(clus);
        if fringe_start > fringe_end
            break;
        end
    end
    IDX(clus) = clusterId;
end
