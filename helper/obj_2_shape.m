function [Shape] = obj_2_shape(objMesh)
% This function takes a wavefront object file and converts it into a Shape
% file used by all the algorithms
% Input argument:
%    objMesh.vertexPoss: the array of vertex positions
%    objMesh.faceVIds: the array of triangle face indices
% Output argument:
%    Shape.vertexPoss: the re-ordered vertex position array
%    Shape.faceVIds: the re-ordered face vertex index array
%    Shape.meshes: the disconnected compoents
Shape.meshes = parts_from_shape(objMesh, 1e-16);
nv = size(objMesh.vertexPoss, 2);
nf = size(objMesh.faceVIds, 2);
Shape.vertexPoss = zeros(3, nv);
Shape.faceVIds = zeros(3, nf);

nv = 0;
nf = 0;
for i = 1:length(Shape.meshes)
    mesh = Shape.meshes{i};
    nv_mesh = length(mesh.vertexIds);
    nf_mesh = size(mesh.faceVIds, 2);
    Shape.vertexPoss(:, (nv+1):(nv+nv_mesh)) =...
        objMesh.vertexPoss(:, mesh.vertexIds);
    Shape.faceVIds(:, (nf+1):(nf+nf_mesh)) =...
        mesh.faceVIds + nv;
    Shape.meshes{i}.vertexIds = (nv+1):(nv+nv_mesh);
    nv = nv + nv_mesh;
    nf = nf + nf_mesh;
end
Shape.has_material = 1;

function [meshes] = parts_from_shape(Shape, err)
%
numVertices = size(Shape.vertexPoss, 2);
ids = zeros(1, numVertices);

knn = min(numVertices, 32);
[TP, DIS] = knnsearch(Shape.vertexPoss', Shape.vertexPoss', 'k', knn);

id = 1;
for i = 1:numVertices
    if ids(i) > 0
        continue;
    end
    tp = TP(i,:);
    dis = DIS(i,:);
    tp = tp(find(dis < err));
    
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
